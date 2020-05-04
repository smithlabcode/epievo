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
#include <cassert>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <random>
#include <functional>
#include <iomanip>

#include "SingleSiteSampler.hpp"
#include "PhyloTreePreorder.hpp"
#include "Path.hpp"
#include "ParamEstimation.hpp"
#include "EndCondSampling.hpp"

#include "epievo_utils.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::string;
using std::pair;
using std::begin;
using std::end;

using std::bind;
using std::placeholders::_1;
using std::multiplies;
using std::accumulate;
using std::runtime_error;
using std::function;
using std::uniform_real_distribution;
using std::exponential_distribution;

static const size_t n_triples = 8;

////////////////////////////////////////////////////////////////////////////////
///////////////////           HELPER CLASSES           /////////////////////////
////////////////////////////////////////////////////////////////////////////////

struct FelsHelper {
  FelsHelper() {}
  FelsHelper(const vector<vector<double> > &_p,
             const vector<double> &_q) : p(_p), q(_q) {}
  vector<vector<double> > p;
  vector<double> q;
};

////////////////////////////////////////////////////////////////////////////////
///////////////////           UPWARD PRUNING           /////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* seg_info gives 2 muation rates: for patterns L0R and L1R. The "q"
 * and "p" values are as defined for Felsenstein's algorithm.
 */
static void
process_interval(const SegmentInfo &si, const vector<double> &q,
                 vector<double> &p) {
  two_by_two P; // transition matrix
  continuous_time_trans_prob_mat(si.rate0, si.rate1, si.len, P);

  // p <- P*q
  p = { P(0, 0)*q[0] + P(0, 1)*q[1], P(1, 0)*q[0] + P(1, 1)*q[1] };
}


static void
process_branch_above(const vector<SegmentInfo> &seg_info, FelsHelper &fh) {

  const size_t n_intervals = seg_info.size();
  fh.p.resize(n_intervals, {0.0, 0.0});

  // directly use q for final interval (back of p[node_id])
  assert(n_intervals > 0);
  process_interval(seg_info.back(), fh.q, fh.p.back());

  // iterate backwards/upwards with: q[i-1] = p[i]
  for (size_t i = n_intervals - 1; i > 0; --i)
    process_interval(seg_info[i-1], fh.p[i], fh.p[i-1]);

  // swap(fh.p, p); // assign the computed p to fh
}


/* first: computes the "q" depending on whether or not node_id
 * indicates a leaf node. If not at leaf node, then the children must
 * be processed, which requires that their "p" values have already
 * been set. So this function must be called in reverse pre-order.
 *
 * next: compute "p" by calling "process_branch_above"
 */
static void
pruning_branch(const TreeHelper &th, const size_t node_id, const Path &path,
               const vector<SegmentInfo> &seg_info, vector<FelsHelper> &fh) {

  /* first calculate q */
  vector<double> q(2, 1.0);
  if (th.is_leaf(node_id)) {
    const bool leaf_state = path.end_state();
    q[0] = (leaf_state == false) ? 1.0 : 0.0;
    q[1] = (leaf_state == true)  ? 1.0 : 0.0;
  }
  else {
    for (ChildSet c(th.subtree_sizes, node_id); c.good(); ++c) {
      q[0] *= fh[*c].p.front()[0];
      q[1] *= fh[*c].p.front()[1];
    }
  }
  fh[node_id].q.swap(q); // assign computed q to the fh

  /* now calculate p if we are not at root */
  if (!th.is_root(node_id))
    process_branch_above(seg_info, fh[node_id]);
}

/* Iterates over all nodes, in reverse pre-order, and for each it sets
 * the "p" and "q" from Felsenstein's algorithm, requiring for each
 * node that the children have already been processed (hence reversing
 * the order).
 */
void
pruning(const TreeHelper &th, const size_t site_id,
        const vector<vector<Path> > &paths,
        const vector<vector<SegmentInfo> > &seg_info, vector<FelsHelper> &fh) {

  // avoid recursion by iterating backwards ensuring all child nodes
  // are processed before their parents
  fh.resize(th.n_nodes);
  for (size_t i = th.n_nodes; i > 0; --i) {
    const size_t node_id = i - 1;
    pruning_branch(th, node_id, paths[site_id][node_id], seg_info[node_id], fh);
  }
}

////////////////////////////////////////////////////////////////////////////////
///////////////            DOWNWARD SAMPLING        ////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* Compute posterior probability of state 0 at root node (i.e. the
 * init_state of any/all children) for particular site using
 * information from the upward stage of Felsenstein's algorithm.
 */
static double
root_post_prob0(const bool left_st, const bool right_st,
                const two_by_two &horiz_tr_prob,
                const vector<double> &q) {

  const double p0 = (horiz_tr_prob(left_st, 0)*horiz_tr_prob(0, right_st))*q[0];
  const double p1 = (horiz_tr_prob(left_st, 1)*horiz_tr_prob(1, right_st))*q[1];

  return p0/(p0 + p1);
}


/* downward sampling on a single branch, given the initial state */
static void
downward_sampling_branch(const vector<SegmentInfo> &seg_info,
                         const FelsHelper &fh,
                         const size_t start_state,
                         const double branch_length,
                         std::mt19937 &gen,
                         uniform_real_distribution<double> &unif,
                         vector<exponential_distribution<double> > &exp,
                         Path &sampled_path, double &log_prob) {

  sampled_path.init_state = start_state;
  sampled_path.tot_time = branch_length;
  sampled_path.jumps.clear();

  size_t prev_state = sampled_path.init_state;
  double time_passed = 0.0;
  const size_t n_intervals = seg_info.size();
  for (size_t i = 0; i < n_intervals; ++i) {

    const TwoStateCTMarkovModel ctmm(seg_info[i].rate0, seg_info[i].rate1);
    const double PT0 = ctmm.get_trans_prob(seg_info[i].len, prev_state, 0);

    const double p0 =
      PT0*((i == n_intervals - 1) ?
           fh.q[0] : fh.p[i + 1][0])/fh.p[i][prev_state];

    const size_t sampled_state = (unif(gen) > p0);
    log_prob += (sampled_state == 0) ? log(p0) : log(1.0 - p0);

    end_cond_sample_forward_rejection(exp[seg_info[i].trip0],
                                      exp[seg_info[i].trip1],
                                      prev_state, sampled_state,
                                      seg_info[i].len, gen, sampled_path.jumps,
                                      time_passed);
    log_prob -= log(ctmm.get_trans_prob(seg_info[i].len, prev_state,
                                        sampled_state));

    // prepare for next interval
    time_passed += seg_info[i].len;
    prev_state = sampled_state;
  }
}


/* Iterates over all nodes in pre-order, and for each branch it samples new
 * paths from top to bottom, requiring starting state has been determined.
 */
void
downward_sampling(const EpiEvoModel &mod, const TreeHelper &th,
                  const bool left_root_state, const bool mid_root_state,
                  const bool right_root_state,
                  const vector<vector<SegmentInfo> > &seg_info,
                  const vector<FelsHelper> &fh,
                  std::mt19937 &gen, uniform_real_distribution<double> &unif,
                  vector<exponential_distribution<double> > &exp,
                  const bool sample_root,
                  vector<Path> &proposed_path, double &log_prob) {

  // compute posterior probability at root node
  const double root_p0 = root_post_prob0(left_root_state, right_root_state,
                                         mod.T, fh[0].q);
  if (sample_root) {
    proposed_path.front().init_state = (unif(gen) > root_p0);
    log_prob = proposed_path.front().init_state ?
      log(1.0 - root_p0) : log(root_p0);
  } else {
    proposed_path.front().init_state = mid_root_state;
  }

  for (size_t node_id = 1; node_id < th.n_nodes; ++node_id) {
    const size_t start_state = proposed_path[th.parent_ids[node_id]].end_state();
    downward_sampling_branch(seg_info[node_id], fh[node_id], start_state,
                             th.branches[node_id], gen, unif, exp,
                             proposed_path[node_id], log_prob);
  }
}


////////////////////////////////////////////////////////////////////////////////
///////////////        Compute acceptance rate         /////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* Compute likelihood of root state at given site and neighboring states */
static double
root_prior_lh(const size_t l, const size_t m, const size_t r,
              const two_by_two &horiz_tr_prob) {

  const double p = horiz_tr_prob(l, m)*horiz_tr_prob(m, r);
  return p;
}

/* Counterpart of downward_sampling_branch */
static double
proposal_prob_branch(const vector<SegmentInfo> &seg_info,
                     const FelsHelper &fh, const Path &path) {

  const size_t n_intervals = seg_info.size();

  size_t start_state = path.init_state, end_state = path.init_state;
  double end_time = 0.0;
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
    const double PT0 = ctmm.get_trans_prob(seg_info[i].len, start_state, 0);
    log_prob -= log(ctmm.get_trans_prob(seg_info[i].len, start_state, end_state));

    // p0 = P_v(j, k) x q_k(v)/p_j(v) [along a branch, q[i]=p[i+1]
    const double p0 = PT0/fh.p[i][start_state]*((i == n_intervals-1) ?
                                                fh.q[0] : fh.p[i+1][0]);
    log_prob += (end_state == 0) ? log(p0) : log(1.0 - p0);
    assert(std::isfinite(log_prob));

    // prepare for next interval
    start_jump = end_jump;
    start_state = end_state;
  }
  return log_prob;
}


/* compute proposal prob with a single path (along all nodes) as input*/
double
proposal_prob(const double(&triplet_rates)[8], const TreeHelper &th,
              const bool rt_left_st, const bool rt_right_st,
              const two_by_two &horiz_trans_prob,
              const vector<FelsHelper> &fh,
              const vector<vector<SegmentInfo> > &seg_info,
              const vector<Path> &the_path, const bool sample_root) {

  // compute posterior probability of state 0 at root node
  double log_prob = 0.0;
  if (sample_root) {
    const double root_p0 =
    root_post_prob0(rt_left_st, rt_right_st, horiz_trans_prob, fh[0].q);
    log_prob += the_path[1].init_state ? log(1.0 - root_p0) : log(root_p0);
  }

  // process the paths above each node (except the root)
  for (size_t node_id = 1; node_id < th.n_nodes; ++node_id) {
    log_prob += proposal_prob_branch(seg_info[node_id], fh[node_id],
                                     the_path[node_id]);
  }
  assert(std::isfinite(log_prob));
  return log_prob;
}


static double
log_likelihood(const double(&rates)[8],
               const vector<double> &J, const vector<double> &D,
               const double(&log_rates)[8]) {
  static const size_t n_triples = 8;
  double r = 0.0;
  for (size_t i = 0; i < n_triples; ++i) {
    r += J[i]*log_rates[i] - D[i]*rates[i];
  }
  return r;
}


/* Currently, emission probabilities are not considered. */
double
path_log_likelihood(const EpiEvoModel &mod, const vector<Path> &l,
                    const vector<Path> &m, const vector<Path> &r,
                    const double(&log_rates)[8]) {
  vector<double> D(n_triples, 0.0), J(n_triples, 0.0);

  /* calculate likelihood involving root states */
  double llh = root_prior_lh(l[1].init_state,
                             m[1].init_state, r[1].init_state, mod.T);

  for (size_t i = 1; i < m.size(); ++i)
    add_sufficient_statistics(l[i], m[i], r[i], J, D);
  llh += log_likelihood(mod.triplet_rates, J, D, log_rates);

  return llh;
}

/* Only when you know J and D have proper size */
static double
path_log_likelihood(const EpiEvoModel &mod, const vector<Path> &l,
                    const vector<Path> &m, const vector<Path> &r,
                    vector<double> &J, vector<double> &D,
                    const double(&log_rates)[8]) {
  fill_n(begin(D), n_triples, 0.0);
  fill_n(begin(J), n_triples, 0.0);

  /* calculate likelihood involving root states */
  double llh =
    root_prior_lh(l[1].init_state, m[1].init_state, r[1].init_state, mod.T);

  for (size_t i = 1; i < m.size(); ++i)
    add_sufficient_statistics(l[i], m[i], r[i], J, D);
  llh += log_likelihood(mod.triplet_rates, J, D, log_rates);

  return llh;
}


/* compute acceptance rate */
/* Currently, emission probabilities are not considered. */
static double
log_accept_rate(const EpiEvoModel &mod, const TreeHelper &th,
                const size_t site_id, const vector<vector<Path> > &paths,
                double &llh_l, double &llh_m, double &llh_r,
                vector<double> &J, vector<double> &D,
                const vector<FelsHelper> &fh,
                const vector<vector<SegmentInfo> > &seg_info,
                const vector<Path> &proposed_path,
                const double update_proposal, const bool sample_root,
                const double(&log_rates)[8]) {

  assert(site_id > 0 && site_id < paths.size());

  const bool rt_left_st = paths[site_id-1][1].init_state;
  const bool rt_right_st = paths[site_id+1][1].init_state;

  const double orig_proposal = proposal_prob(mod.triplet_rates, th,
                                             rt_left_st, rt_right_st,
                                             mod.T, fh, seg_info,
                                             paths[site_id], sample_root);

  double llr = orig_proposal - update_proposal;

  const double llh_l_orig = llh_l;
  const double llh_r_orig = llh_r;

  if (site_id > 1)
    llh_l = path_log_likelihood(mod, paths[site_id-2], paths[site_id-1],
                                proposed_path, J, D, log_rates);
  llh_m = path_log_likelihood(mod, paths[site_id-1], proposed_path,
                              paths[site_id+1], J, D, log_rates);
  if (site_id < paths.size() - 2)
    llh_r = path_log_likelihood(mod, proposed_path, paths[site_id+1],
                                paths[site_id+2], J, D, log_rates);

  llr += (llh_l + llh_r - llh_l_orig - llh_r_orig);
  return llr;
}

////////////////////////////////////////////////////////////////////////////////
///////////////               MCMC setup               /////////////////////////
////////////////////////////////////////////////////////////////////////////////

SingleSiteSampler::SingleSiteSampler(const size_t n_burn_in,
                                     const size_t n_batch) {
  SAMPLE_ROOT = false;
  burn_in = n_burn_in;
  batch = n_batch;

  unif = uniform_real_distribution<double> (0.0, 1.0);
  exp_distr.resize(8);
}

void
SingleSiteSampler::reset(const EpiEvoModel &the_model,
                         vector<vector<Path> > &paths) {
  const size_t n_sites = paths.size();
  const size_t n_nodes = paths.front().size();

  J.resize(n_triples, 0.0);
  D.resize(n_triples, 0.0);
  proposed_path.resize(n_nodes);
  tri_llh.resize(n_sites, 0.0);

  // update exponential distributions
  for(size_t i = 0; i < the_model.n_triplets; i++)
    exp_distr[i] = exponential_distribution<double>(the_model.triplet_rates[i]);

  // pre-compute log(rates)
  std::transform(begin(the_model.triplet_rates),
                 end(the_model.triplet_rates),
                 begin(log_rates),
                 static_cast<double(*)(double)>(log));

  for (size_t site_id = 1; site_id < n_sites - 1; ++site_id) {
    tri_llh[site_id] = path_log_likelihood(the_model, paths[site_id-1],
                                           paths[site_id], paths[site_id+1],
                                           log_rates);
  }
}

////////////////////////////////////////////////////////////////////////////////
///////////////          single MCMC iteration         /////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* Metropolis-Hastings sampling at single site */
bool
SingleSiteSampler::Metropolis_Hastings_site(const EpiEvoModel &the_model,
                                            const TreeHelper &th,
                                            const size_t site_id,
                                            vector<vector<Path> > &paths,
                                            std::mt19937 &gen) {
  // get rates and lengths each interval [seg_info: node x segs]
  vector<vector<SegmentInfo> > seg_info(th.n_nodes);
  for (size_t node_id = 1; node_id < th.n_nodes; ++node_id) {
    collect_segment_info(the_model.triplet_rates,
                         paths[site_id - 1][node_id],
                         paths[site_id + 1][node_id], seg_info[node_id]);
  }

  // upward pruning and downward sampling [fh: one for each node]
  vector<FelsHelper> fh;
  pruning(th, site_id, paths, seg_info, fh);

  double proposal_log_prob = 0.0;
  downward_sampling(the_model, th,
                    paths[site_id-1][1].init_state,
                    paths[site_id][1].init_state,
                    paths[site_id+1][1].init_state,
                    seg_info, fh, gen, unif, exp_distr,
                    SAMPLE_ROOT,
                    proposed_path, proposal_log_prob);

  // acceptance rate
  double llh_l_prop = tri_llh[site_id - 1];
  double llh_m_prop = tri_llh[site_id];
  double llh_r_prop = tri_llh[site_id + 1];

  const double log_acc_rate =
    log_accept_rate(the_model, th, site_id, paths,
                    llh_l_prop, llh_m_prop, llh_r_prop, J, D,
                    fh, seg_info, proposed_path, proposal_log_prob,
                    SAMPLE_ROOT, log_rates);

  std::uniform_real_distribution<double> unif(0.0, 1.0);
  const double u = unif(gen);

  bool accepted = false;
  if (log_acc_rate >= 0 || u < exp(log_acc_rate))
    accepted = true;

  // if accepted, replace old path with proposed one.
  if (accepted) {
    std::swap(paths[site_id], proposed_path);
    tri_llh[site_id - 1] = llh_l_prop;
    tri_llh[site_id] = llh_m_prop;
    tri_llh[site_id + 1] = llh_r_prop;
  }

  return accepted;
}

size_t
SingleSiteSampler::single_iteration(const EpiEvoModel &the_model,
                                    const TreeHelper &th,
                                    vector<vector<Path> > &paths,
                                    std::mt19937 &gen) {

  size_t n_accepted = 0;
  for (size_t site_id = 1; site_id < paths.size() - 1; ++site_id)
    n_accepted += Metropolis_Hastings_site(the_model, th, site_id, paths, gen);
  return n_accepted;
}

void
SingleSiteSampler::run_mcmc(const EpiEvoModel &the_model,
                            const TreeHelper &th,
                            vector<vector<Path> > &paths,
                            std::mt19937 &gen,
                            vector<vector<double> > &J_all_sites,
                            vector<vector<double> > &D_all_sites,
                            double &acc_rate) {

  // Burning
  for(size_t burnin_itr = 0; burnin_itr < burn_in; burnin_itr++)
    single_iteration(the_model, th, paths, gen);

  // Prepare for taking MCMC samples
  J_all_sites.resize(th.n_nodes);
  D_all_sites.resize(th.n_nodes);
  for(size_t b = 1; b < th.n_nodes; ++b) {
    J_all_sites[b].clear();
    J_all_sites[b].resize(n_triples, 0.0);
    D_all_sites[b].clear();
    D_all_sites[b].resize(n_triples, 0.0);
  }

  size_t n_accepted = 0;
  vector<vector<double> > J_one_site, D_one_site;

  for (size_t mcmc_itr = 0; mcmc_itr < batch; mcmc_itr++) {
    n_accepted += single_iteration(the_model, th, paths, gen);

    /* CALCULATE SUFFICIENT STATS */
    get_sufficient_statistics(paths, J_one_site, D_one_site);

    for (size_t b = 1; b < th.n_nodes; b++)
      for (size_t i = 0; i < n_triples; i++) {
        J_all_sites[b][i] += J_one_site[b][i];
        D_all_sites[b][i] += D_one_site[b][i];
      }
    }

  /* CALCULATE BATCH AVERAGE */
  for (size_t b = 1; b < th.n_nodes; ++b)
    for (size_t i = 0; i < n_triples; i++) {
      J_all_sites[b][i] /= batch;
      D_all_sites[b][i] /= batch;
    }

  /* CALCULATE ACCEPTANCE RATE */
  acc_rate = static_cast<double>(n_accepted)/(batch*(paths.size() - 2));
}
