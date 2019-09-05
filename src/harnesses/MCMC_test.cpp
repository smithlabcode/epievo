/* Copyright (C) 2019 University of Southern California
 *                    Xiaojing Ji, Jianghan Qu and Andrew D Smith
 *
 * Author: Andrew D. Smith, Jianghan Qu and Xiaojing Ji
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA
 */

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <random>
#include <functional>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "Path.hpp"
#include "EpiEvoModel.hpp"
#include "EndCondSampling.hpp"
#include "ContinuousTimeMarkovModel.hpp"
#include "TripletSampler.hpp"
#include "GlobalJump.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::runtime_error;
using std::numeric_limits;
using std::begin;
using std::end;

using std::exponential_distribution;
using std::placeholders::_1;
using std::function;


template <typename T>
using disc_distr = std::discrete_distribution<T>;

/* below: calls functions declared in Path.hpp */
void
add_sufficient_statistics(const vector<Path> &paths,
                          vector<double> &J, vector<double> &D) {
  // iterate over sites with valid triples (i.e. not the first and last)
  const size_t n_sites = paths.size();
  for (size_t i = 1; i < n_sites - 1; ++i)
    add_sufficient_statistics(paths[i-1], paths[i], paths[i+1], J, D);
}

struct SegmentInfo {
  SegmentInfo() {}
  SegmentInfo(const double r0, const double r1, const double l) :
    rate0(r0), rate1(r1), len(l) {}
  double rate0;
  double rate1;
  double len;
};

struct FelsHelper {
  FelsHelper() {}
  FelsHelper(const std::vector<std::vector<double> > &_p,
             const std::vector<double> &_q) : p(_p), q(_q) {}
  std::vector<std::vector<double> > p;
  std::vector<double> q;
};

struct Environment {
  vector<bool> left; // states on the left
  vector<bool> right; // states on the right
  vector<double> breaks; // times for state change, including tot_time
  double tot_time;

  Environment(const Path &pl, const Path &pr) {
    bool sl = pl.init_state;
    bool sr = pr.init_state;
    size_t i = 0;
    size_t j = 0;
    tot_time = pl.tot_time;
    while (i < pl.jumps.size() && j < pr.jumps.size()) {
      left.push_back(sl);
      right.push_back(sr);
      if (pl.jumps[i] < pr.jumps[j]) {
        breaks.push_back(pl.jumps[i]);
        ++i;
        sl = !sl;
      }
      else if (pr.jumps[j] < pl.jumps[i]) {
        breaks.push_back(pr.jumps[j]);
        ++j;
        sr = !sr;
      }
      else {
        breaks.push_back(pr.jumps[j]);
        ++j;
        ++i;
        sl = !sl;
        sr = !sr;
      }
    }
    while (i < pl.jumps.size()) {
      left.push_back(sl);
      right.push_back(sr);
      breaks.push_back(pl.jumps[i]);
      ++i;
      sl = !sl;
    }
    while (j < pr.jumps.size()) {
      left.push_back(sl);
      right.push_back(sr);
      breaks.push_back(pr.jumps[j]);
      ++j;
      sr = !sr;
    }
    if (breaks.empty() || breaks.back() < tot_time) {
      left.push_back(sl);
      right.push_back(sr);
      breaks.push_back(tot_time);
    }
  }
};


inline size_t
triple2idx(const bool i, const bool j, const bool k) {
  return i*4 + j*2 + k;
}


/* collect rates and interval lengths */
static void
collect_segment_info(const vector<double> &triplet_rates,
                     const Path &l, const Path &r,
                     vector<SegmentInfo> &seg_info) {
  Environment env(l, r);
  const size_t n_intervals = env.left.size();
  seg_info = vector<SegmentInfo>(n_intervals);

  for (size_t i = 0; i < n_intervals; ++i) {
    const size_t pattern0 = triple2idx(env.left[i], false, env.right[i]);
    const size_t pattern1 = triple2idx(env.left[i], true, env.right[i]);
    seg_info[i] = SegmentInfo(triplet_rates[pattern0], triplet_rates[pattern1],
                              env.breaks[i] - (i == 0 ? 0.0 : env.breaks[i-1]));
    assert(seg_info[i].len > 0.0);
  }
}


////////////////////////////////////////////////////////////////////////////////
///////////////////           UPWARD PRUNING           /////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* seg_info gives 2 muation rates: for patterns L0R and L1R. The "q"
 * and "p" values are as defined for Felsenstein's algorithm.
 */
static void
process_interval(const SegmentInfo &si, const vector<double> &q,
                 vector<double> &p) {
  vector<vector<double> > P; // transition matrix
  continuous_time_trans_prob_mat(si.rate0, si.rate1, si.len, P);

  // p <- P*q
  p = { P[0][0]*q[0] + P[0][1]*q[1], P[1][0]*q[0] + P[1][1]*q[1] };
}


static void
process_above(const vector<SegmentInfo> &seg_info, FelsHelper &fh) {

  const size_t n_intervals = seg_info.size();
  vector<vector<double> > p(n_intervals, {0.0, 0.0});

  // directly use q for final interval (back of p[node_id])
  assert(n_intervals > 0);
  process_interval(seg_info.back(), fh.q, p.back());

  // iterate backwards/upwards with: q[i-1] = p[i]
  for (size_t i = n_intervals - 1; i > 0; --i)
    process_interval(seg_info[i-1], p[i], p[i-1]);

  swap(fh.p, p); // assign the computed p to fh
}


/* first: computes the "q" depending on whether or not node_id
 * indicates a leaf node. If not at leaf node, then the children must
 * be processed, which requires that their "p" values have already
 * been set. So this function must be called in reverse pre-order.
 *
 * next: compute "p" by calling "process_above"
 */
static void
pruning(const Path &the_path, const vector<SegmentInfo> &seg_info,
        FelsHelper &fh) {

  /* first calculate q */
  vector<double> q = { 1.0, 1.0 };
  const bool leaf_state = the_path.end_state();
  q[0] = (leaf_state == false) ? 1.0 : 0.0;
  q[1] = (leaf_state == true)  ? 1.0 : 0.0;
  fh.q.swap(q); // assign computed q to the fh

  process_above(seg_info, fh);
}

////////////////////////////////////////////////////////////////////////////////
///////////////            DOWNWARD SAMPLING        ////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* Compute posterior probability of state 0 at root node (i.e. the
 * init_state of any/all children) for particular site using
 * information from the upward stage of Felsenstein's algorithm.
 */
static double
root_post_prob0(const size_t site_id, const vector<Path> &the_paths,
                const vector<vector<double> > &horiz_tr_prob,
                const vector<double> &q) {

  const size_t left_st = the_paths[site_id - 1].init_state;
  const size_t right_st = the_paths[site_id + 1].init_state;

  const double p0 = (horiz_tr_prob[left_st][0]*horiz_tr_prob[0][right_st])*q[0];
  const double p1 = (horiz_tr_prob[left_st][1]*horiz_tr_prob[1][right_st])*q[1];

  return p0/(p0 + p1);
}


/* downward sampling on a single branch, given the initial state */
static void
downward_sampling(const vector<SegmentInfo> &seg_info,
                  const FelsHelper &fh,
                  const bool start_state,
                  const double branch_length,
                  std::mt19937 &gen,
                  Path &sampled_path) {

  sampled_path.init_state = start_state;
  sampled_path.tot_time = branch_length;
  sampled_path.jumps.clear();

  std::uniform_real_distribution<double> unif(0.0, 1.0);

  bool prev_state = sampled_path.init_state;
  double time_passed = 0.0;
  for (size_t i = 0; i < seg_info.size(); ++i) {

    const TwoStateCTMarkovModel ctmm(seg_info[i].rate0, seg_info[i].rate1);
    const double PT0 = ctmm.get_trans_prob(seg_info[i].len, prev_state, 0);

    const double p0 =
      PT0*((i == seg_info.size() - 1) ? fh.q[0] : fh.p[i + 1][0])/fh.p[i][prev_state];

    const bool sampled_state = (unif(gen) > p0);

    end_cond_sample_forward_rejection(ctmm, prev_state, sampled_state,
                                      seg_info[i].len, gen, sampled_path.jumps,
                                      time_passed);

    // prepare for next interval
    time_passed += seg_info[i].len;
    prev_state = sampled_state;
  }
}


////////////////////////////////////////////////////////////////////////////////
///////////////        Compute acceptance rate         /////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* Compute likelihood of root state at given site and neighboring states */
static double
root_prior_lh(const size_t l, const size_t m, const size_t r,
              const vector<vector<double> > &horiz_tr_prob) {
  const double p = horiz_tr_prob[l][m]*horiz_tr_prob[m][r];
  return p;
}

inline bool
complement_state(const bool x) {
  return !x;
}

/* Counterpart of downward_sampling_branch */
static double
proposal_prob(const vector<SegmentInfo> &seg_info,
              const FelsHelper &fh, const Path &path) {

  const size_t n_intervals = seg_info.size();

  bool start_state = path.init_state, end_state = path.init_state;
  double start_time = 0.0, end_time = 0.0;
  size_t start_jump = 0, end_jump = 0;

  double log_prob = 0.0;
  for (size_t i = 0; i < n_intervals; ++i) {

    // get the end_time by adding the segment length to the start_time
    end_time += seg_info[i].len;

    // get the end_jump as the first jump occurring after end_time
    while (end_jump < path.jumps.size() && path.jumps[end_jump] < end_time)
      ++end_jump;

    // get end_state based on number of jumps between start_time and end_time
    if ((end_jump - start_jump) % 2 == 1)
      end_state = complement_state(end_state);

    // calculate the probability for the end-conditioned path
    const TwoStateCTMarkovModel ctmm(seg_info[i].rate0, seg_info[i].rate1);

    const double interval_prob =
        end_cond_sample_prob(ctmm, path.jumps, start_state, end_state,
                             start_time, end_time, start_jump, end_jump);
    log_prob += interval_prob;

    const double PT0 = ctmm.get_trans_prob(seg_info[i].len, start_state, 0);

    // p0 = P_v(j, k) x q_k(v)/p_j(v) [along a branch, q[i]=p[i+1]
    const double p0 = PT0/fh.p[i][start_state]*((i == n_intervals-1) ?
                                                fh.q[0] : fh.p[i+1][0]);

    log_prob += (end_state == 0) ? log(p0) : log(1.0 - p0);
    assert(std::isfinite(log_prob));

    // prepare for next interval
    start_jump = end_jump;
    start_time = end_time;
    start_state = end_state;
  }
  return log_prob;
}


/* compute proposal prob */
static double
proposal_prob(const vector<double> &triplet_rates,
              const double evo_time,
              const size_t site_id,
              const vector<Path> &paths,
              const vector<vector<double> > &horiz_trans_prob,
              const FelsHelper &fh,
              const vector<SegmentInfo> &seg_info,
              const Path &the_path) {

  // compute posterior probability of state 0 at root node
  const double root_p0 =
    root_post_prob0(site_id, paths, horiz_trans_prob, fh.p.front());
  const double log_prob =
    (the_path.init_state ? log(1.0 - root_p0) : log(root_p0)) +
    proposal_prob(seg_info, fh, the_path);

  return log_prob;
}


static double
log_likelihood(const vector<double> &rates,
               const vector<double> &J, const vector<double> &D) {
  static const size_t n_triples = 8;
  double r = 0.0;
  for (size_t i = 0; i < n_triples; ++i) {
    r += J[i]*log(rates[i]) - D[i]*rates[i];
  }
  return r;
}


/* compute acceptance rate */
double
log_accept_rate(const EpiEvoModel &mod, const double evo_time,
                const size_t site_id,
                const vector<Path> &paths,
                const FelsHelper &fh,
                const vector<SegmentInfo> &seg_info,
                const Path &proposed_path) {

  static const size_t n_triples = 8;

  // ADS: this is unfortunate; we need to slice/transpose the original
  // paths because of how they are organized
  Path original(paths[site_id]);
  const double orig_proposal =
    proposal_prob(mod.triplet_rates, evo_time, site_id, paths,
                  mod.init_T, fh, seg_info, original);

  const double update_proposal =
    proposal_prob(mod.triplet_rates, evo_time, site_id, paths,
                  mod.init_T, fh, seg_info, proposed_path);

  assert(site_id > 0);

  double llr = orig_proposal - update_proposal;

  /* calculate likelihood involving root states */
  const size_t rt_l = paths[site_id-1].init_state;
  const size_t rt_r = paths[site_id+1].init_state;
  const size_t rt_orig = paths[site_id].init_state;
  const size_t rt_prop = proposed_path.init_state;
  llr += (log(root_prior_lh(rt_l, rt_prop, rt_r, mod.init_T)) -
          log(root_prior_lh(rt_l, rt_orig, rt_r, mod.init_T)));

  /* calculate likelihood involving internal intervals */
  vector<double> D_orig(n_triples), J_orig(n_triples);
  vector<double> D_prop(n_triples), J_prop(n_triples);

  // sufficient stats for current pentet using original path at mid
  vector<Path>::const_iterator opth(begin(paths) + site_id);
  fill_n(begin(D_orig), n_triples, 0.0);
  fill_n(begin(J_orig), n_triples, 0.0);

  if (site_id > 1)
    add_sufficient_statistics(*(opth - 2), *(opth - 1), *opth, J_orig, D_orig);
  add_sufficient_statistics(*(opth - 1), *opth, *(opth + 1), J_orig, D_orig);
  if (site_id < paths.size() - 2)
    add_sufficient_statistics(*opth, *(opth + 1), *(opth + 2), J_orig, D_orig);

  // sufficient stats for current pentet using proposed path at mid
  fill_n(begin(D_prop), n_triples, 0.0);
  fill_n(begin(J_prop), n_triples, 0.0);

  if (site_id > 1)
    add_sufficient_statistics(*(opth - 2), *(opth - 1), proposed_path,
                              J_prop, D_prop);
  add_sufficient_statistics(*(opth - 1), proposed_path, *(opth + 1),
                            J_prop, D_prop);
  if (site_id < paths.size() - 2)
    add_sufficient_statistics(proposed_path, *(opth + 1), *(opth + 2),
                              J_prop, D_prop);

  // add difference in log-likelihood for the proposed vs. original
  // to the Hastings ratio
  llr += (log_likelihood(mod.triplet_rates, J_prop, D_prop) -
          log_likelihood(mod.triplet_rates, J_orig, D_orig));
  return llr;
}

////////////////////////////////////////////////////////////////////////////////
///////////////          single MCMC iteration         /////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* Metropolis-Hastings sampling at single site */
bool
Metropolis_Hastings_site(const EpiEvoModel &the_model, const double &evo_time,
                         const size_t site_id, vector<Path> &paths,
                         std::mt19937 &gen, Path &proposed_path) {

  // get rates and lengths each interval
  vector<SegmentInfo> seg_info;
  collect_segment_info(the_model.triplet_rates,
                       paths[site_id - 1], paths[site_id + 1], seg_info);

  // upward pruning and downward sampling
  FelsHelper fh;
  pruning(paths[site_id], seg_info, fh);
  downward_sampling(seg_info, fh, paths[site_id].init_state,
                    evo_time, gen, proposed_path);

  // acceptance rate
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  const double u = unif(gen);
  const double log_acc_rate =
    log_accept_rate(the_model, evo_time, site_id, paths,
                    fh, seg_info, proposed_path);

  bool accepted = false;
  if (log_acc_rate >= 0 || u < exp(log_acc_rate))
    accepted = true;

  // if accepted, replace old path with proposed one.
  if (accepted)
    paths[site_id] = proposed_path;

  return accepted;
}

static void
collect_init_sequences(const vector<Path> &paths, vector<bool> &seq) {
  seq = vector<bool>(paths.size());
  for (size_t i = 0; i < paths.size(); ++i)
    seq[i] = paths[i].init_state;
}

static void
collect_end_sequences(const vector<Path> &paths, vector<bool> &seq) {
  seq = vector<bool>(paths.size());
  for (size_t i = 0; i < paths.size(); i++)
    seq[i] = paths[i].end_state();
}

/* This function does the sampling for an individual change in the
 * state sequence
 */
static void
sample_jump(const EpiEvoModel &the_model, const double total_time,
            std::mt19937 &gen, TripletSampler &ts,
            vector<GlobalJump> &the_path, double &time_value) {

  static const size_t n_triplets = 8;

  // triplet_count = c_{ijk} for current sequence (encoded in the
  // TripletSampler object)
  vector<size_t> triplet_counts;
  ts.get_triplet_counts(triplet_counts);

  // holding_rate = c_{ijk}*lambda_{ijk}
  const double holding_rate =
    std::inner_product(begin(triplet_counts), end(triplet_counts),
                       begin(the_model.triplet_rates), 0.0);
  // sample a holding time = time until next state change
  std::exponential_distribution<double> exp_distr(holding_rate);
  const double holding_time = std::max(exp_distr(gen),
                                       numeric_limits<double>::min());

  // update the current time_value
  time_value += holding_time;

  // if the holding time ends before the total time interval, we can
  // make a change to the state sequence
  if (time_value < total_time) {

    /* first: get a probability distribution for the triplet to change */
    vector<double> triplet_prob(n_triplets, 0.0);
    for (size_t i = 0; i < n_triplets; ++i)
      triplet_prob[i] =
        triplet_counts[i]*the_model.triplet_rates[i]/holding_rate;

    /* next: use that distribution to sample which triplet type at
       which the change will happen */
    disc_distr<size_t> multinom(begin(triplet_prob), end(triplet_prob));
    const size_t context = multinom(gen);

    /* sample a change position having the relevant triplet; this
       changes the TripletSampler data structure to reflect a changed
       state at the position sampled */
    const size_t change_position = ts.random_mutate(context, gen);

    /* add the changed position and change time to the path */
    the_path.push_back(GlobalJump(time_value, change_position));
  }
}


static bool
check_leaf_seq(const vector<bool> &s, const vector<Path> &by_site) {
  bool leaf_is_good = true;
  for (size_t i = 0; i < by_site.size() && leaf_is_good; ++i)
    leaf_is_good = (by_site[i].end_state() == s[i]);
  return leaf_is_good;
}


static void
global_to_local(const vector<bool> &root_seq, const double evo_time,
                const vector<GlobalJump> &gp, vector<Path> &paths) {

  const size_t n_sites = root_seq.size();

  paths = vector<Path>(n_sites);
  for (size_t i = 0; i < n_sites; ++i) {
    paths[i].init_state = root_seq[i];
    paths[i].tot_time = evo_time;
  }

  for (size_t i = 0; i < gp.size(); ++i)
    paths[gp[i].position].jumps.push_back(gp[i].timepoint);
}


static void
forward_simulation(const EpiEvoModel &the_model,
                   const double &curr_branch_len,
                   const size_t n_sites, vector<Path> &paths,
                   std::mt19937 &gen) {

  // sample a root
  vector<bool> root_seq;
  the_model.sample_state_sequence_init(n_sites, gen, root_seq);
  TripletSampler ts(root_seq);

  // generate the global path
  vector<GlobalJump> global_path;
  double time_value = 0;
  while (time_value < curr_branch_len)
    sample_jump(the_model, curr_branch_len, gen, ts, global_path, time_value);

  // convert global to site-specific histories
  global_to_local(root_seq, curr_branch_len, global_path, paths);
}


static bool
end_cond_forward_simulation(const EpiEvoModel &the_model, const double &evo_time,
                            const vector<bool> &root_seq,
                            const vector<bool> &leaf_seq,
                            vector<Path> &paths, std::mt19937 &gen) {

  // start new branch
  const double curr_branch_len = evo_time;
  TripletSampler ts(root_seq);

  vector<GlobalJump> global_path;
  double time_value = 0.0;
  while (time_value < curr_branch_len)
    sample_jump(the_model, curr_branch_len, gen, ts, global_path, time_value);

  vector<bool> s;
  ts.get_sequence(s);

  // Convert global jumps to local paths
  global_to_local(root_seq, curr_branch_len, global_path, paths);

  // leaf sequence agree?
  const bool leaf_is_good = check_leaf_seq(leaf_seq, paths);

  return leaf_is_good;
}


/* generate initial paths by heuristics */
static void
initialize_paths(const vector<bool> &root_seq, const vector<bool> &leaf_seq,
                 const double evo_time, vector<Path> &paths,
                 std::mt19937 &gen) {
  
  paths.resize(root_seq.size());
  
  auto unif =
  bind(std::uniform_real_distribution<double>(0.0, 1.0), std::ref(gen));
  
  for (size_t site_id = 0; site_id < root_seq.size(); ++site_id) {
    paths[site_id].init_state = root_seq[site_id];
    paths[site_id].tot_time = evo_time;
    
    if (root_seq[site_id] != leaf_seq[site_id])
      paths[site_id].jumps.push_back(unif() * evo_time);
  }
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

static void
write_statistics_header(const string outfile) {
  std::ofstream out(outfile.c_str());
  out << "ITR\tSAMPLE\tJ_000\tJ_001\tJ_010\tJ_011\tJ_100\tJ_101\tJ_110\tJ_111\t"
      << "D_000\tD_001\tD_010\tD_011\tD_100\tD_101\tD_110\tD_111" << endl;
}


static void
write_statistics(const string outfile, const size_t itr, const size_t sample,
                 const vector<double> &J, const vector<double> &D) {
  std::ofstream out(outfile, std::ofstream::app);
  out << itr << "\t" << sample << "\t";
  copy(begin(J), begin(J) + 8, std::ostream_iterator<double>(out, "\t"));
  copy(begin(D), begin(D) + 8, std::ostream_iterator<double>(out, "\t"));
  out << endl;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char **argv) {
  try {
    bool VERBOSE = false;

    string outfile;
    string tree_file;

    size_t n_sites = 5;
    size_t batch = 100;
    size_t n_mcmc_batches = 100;

    size_t rng_seed = std::numeric_limits<size_t>::max();

    double evolutionary_time = numeric_limits<double>::lowest();

    ///////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]),
                           "test MCMC against forward-simulation summary stats",
                           "<param>");
    opt_parse.add_opt("n_sites", 'n', "number of sites", false, n_sites);
    opt_parse.add_opt("mcmc_itr", 'i', "number of MCMC iterations",
                      false, n_mcmc_batches);
    opt_parse.add_opt("batch", 'B', "batch size",
                      false, batch);
    opt_parse.add_opt("tree", 't', "Newick format tree file", false, tree_file);
    opt_parse.add_opt("evo-time", 'T', "evolutionary time", true,
                      evolutionary_time);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);;
    opt_parse.add_opt("outfile", 'o', "outfile (prefix)",
                      false, outfile);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() < 1) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    const string param_file(leftover_args[0]);
    ///////////////////////////////////////////////////////////////////////////

    /* (1) READING PARAMETERS FROM FILE */
    if (VERBOSE)
      cerr << "[READING PARAMETERS: " << param_file << "]" << endl;
    EpiEvoModel the_model;
    read_model(param_file, the_model);
    the_model.scale_triplet_rates();
    if (VERBOSE)
      cerr << the_model << endl;

    /* (3) INITIALIZING THE RANDOM NUMBER GENERATOR */
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    std::mt19937 gen(rng_seed);
    if (VERBOSE)
      cerr << "[INITIALIZED RNG (seed=" << rng_seed << ")]" << endl;


    /* (4) GENERATE ONE PATH */
    if (VERBOSE)
      cerr << "[GENERATE A SAMPLE PATH]"<< endl;
    vector<Path> paths; // along multiple branches
    forward_simulation(the_model, evolutionary_time, n_sites, paths, gen);
    // collect sequences
    vector<bool> root_seq;
    collect_init_sequences(paths, root_seq);
    
    vector<bool> leaf_seq;
    collect_end_sequences(paths, leaf_seq);
    if (VERBOSE) {
      cout << "ROOT SEQ: ";
      copy(begin(root_seq), end(root_seq),
           std::ostream_iterator<int>(cout, ""));
      cout << "\nEND SEQ: ";
      copy(begin(leaf_seq), end(leaf_seq),
           std::ostream_iterator<int>(cout, ""));
      cout << endl;
    }
    
    /* (5) FORWARD SIMULATION */
    if (VERBOSE)
      cerr << "[FORWARD SIMULATION]" << endl;
    // forward-simulation output files
    const string fstat = outfile + ".forward";
    write_statistics_header(fstat);

    // sampling
    cerr << "[SAMPLING]" << endl;
    size_t n_forward_samples_collected = 0;
    while (n_forward_samples_collected < batch) {

      vector<Path> sampled_paths;
      bool success = end_cond_forward_simulation(the_model, evolutionary_time,
                                                 root_seq, leaf_seq,
                                                 sampled_paths, gen);
      if (success) {
        vector<double> J(8), D(8);
        add_sufficient_statistics(sampled_paths, J, D);
        write_statistics(fstat, 0, n_forward_samples_collected, J, D);
        n_forward_samples_collected++;
      }
    }

    if (VERBOSE)
      cerr << "[MCMC USING SINGLESITESAMPLER]" << endl;

    // mcmc output files
    const string mcmc_stat = outfile + ".mcmc";
    write_statistics_header(mcmc_stat);

    // distort/randomize paths
    //vector<Path> mcmc_paths(sampled_paths);
    vector<Path> mcmc_paths;
    initialize_paths(root_seq, leaf_seq, evolutionary_time, mcmc_paths, gen);
    vector<Path> proposed_path(n_sites);

    for (size_t batch_id = 0; batch_id < n_mcmc_batches; batch_id++) {
      for (size_t sample_id = 0; sample_id < batch; sample_id++) {
        for (size_t site_id = 1; site_id < n_sites - 1; ++site_id) {
          Metropolis_Hastings_site(the_model, evolutionary_time, site_id,
                                   mcmc_paths, gen, proposed_path[site_id]);
          assert(mcmc_paths[site_id].end_state() == leaf_seq[site_id]);
        }

        // write stats
        vector<double> J(8), D(8);
        add_sufficient_statistics(mcmc_paths, J, D);
        write_statistics(mcmc_stat, batch_id, sample_id, J, D);
      }
    }

  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
