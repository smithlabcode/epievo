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

#include "StateSeq.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::string;
using std::pair;
using std::make_pair;

using std::bind;
using std::placeholders::_1;
using std::multiplies;
using std::runtime_error;

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
process_branch_above(const vector<SegmentInfo> &seg_info, FelsHelper &fh) {

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
 * next: compute "p" by calling "process_branch_above"
 */
static void
pruning_branch(const TreeHelper &th, const size_t site_id, const size_t node_id,
               const vector<Path> &paths,
               const vector<SegmentInfo> &seg_info,
               vector<FelsHelper> &fh) {

  assert(fh.size() == th.n_nodes);

  /* first calculate q */
  vector<double> q(2, 1.0);
  if (th.is_leaf(node_id)) {
    const bool leaf_state = paths[site_id].end_state();
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
        const vector<vector<SegmentInfo> > &seg_info,
        vector<FelsHelper> &fh) {

  // avoid recursion by iterating backwards ensuring all child nodes
  // are processed before their parents
  fh.resize(th.n_nodes);
  for (size_t i = th.n_nodes; i > 0; --i) {
    const size_t node_id = i - 1;
    pruning_branch(th, site_id, node_id, paths[node_id], seg_info[node_id], fh);
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
downward_sampling_branch(const vector<SegmentInfo> &seg_info,
                         const FelsHelper &fh,
                         const size_t start_state,
                         const double branch_length,
                         std::mt19937 &gen,
                         Path &sampled_path, const size_t proposal) {

  sampled_path.init_state = start_state;
  sampled_path.tot_time = branch_length;

  std::uniform_real_distribution<double> unif(0.0, 1.0);

  size_t prev_state = sampled_path.init_state;
  double time_passed = 0;
  const size_t n_intervals = seg_info.size();
  for (size_t i = 0; i < n_intervals; ++i) {

    vector<vector<double> > P;
    continuous_time_trans_prob_mat(seg_info[i].rate0,
                                   seg_info[i].rate1, seg_info[i].len, P);

    const double p0 =
      P[prev_state][0]*((i == n_intervals - 1) ?
                        fh.q[0] : fh.p[i + 1][0])/fh.p[i][prev_state];

    const size_t sampled_state = (unif(gen) > p0);

    const CTMarkovModel ctmm(seg_info[i].rate0, seg_info[i].rate1);

    /* proposal: 0 - direct, 1 - unif, 2 - poisson, 3 - forward */
    if (proposal == 1) { // Uniformization
      end_cond_sample_unif(ctmm, prev_state, sampled_state, seg_info[i].len,
                           gen, sampled_path.jumps, sampled_path.mjumps,
                           time_passed);
      sampled_path.mjump_touched = true;
    } else if (proposal == 2) { // Poisson
      end_cond_sample_Poisson(ctmm, prev_state, sampled_state,
                              seg_info[i].len, gen, sampled_path.jumps,
                              time_passed);
    } else if (proposal == 3) { // Forward
      assert(end_cond_sample_forward_rejection(10000000,
                                               ctmm,
                                               prev_state,
                                               sampled_state,
                                               seg_info[i].len,
                                               gen,
                                               sampled_path.jumps,
                                               time_passed));
    } else {
      end_cond_sample_direct(ctmm, prev_state, sampled_state, seg_info[i].len,
                             gen, sampled_path.jumps, time_passed);
    }

    // prepare for next interval
    time_passed += seg_info[i].len;
    prev_state = sampled_state;
  }
}


/* Iterates over all nodes in pre-order, and for each branch it samples new
 * paths from top to bottom, requiring starting state has been determined.
 */
void
downward_sampling(const TreeHelper &th,
                  const size_t site_id,
                  const vector<vector<Path> > &paths,
                  const vector<vector<double> > &horiz_trans_prob,
                  const vector<vector<SegmentInfo> > &seg_info,
                  const vector<FelsHelper> &fh,
                  std::mt19937 &gen,
                  vector<Path> &proposed_path, const size_t proposal,
                  const bool FIX_ROOT) {

  // compute posterior probability at root node
  const double root_p0 =
    root_post_prob0(site_id, paths[1], horiz_trans_prob, fh[0].q);
  proposed_path = vector<Path>(th.n_nodes);
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  
  proposed_path.front() = Path(paths[1][site_id].init_state, 0.0); // copy root
  if (!FIX_ROOT)
    proposed_path.front() = Path(unif(gen) > root_p0, 0.0); // update root

  for (size_t node_id = 1; node_id < th.n_nodes; ++node_id) {
    const size_t start_state = proposed_path[th.parent_ids[node_id]].end_state();
    downward_sampling_branch(seg_info[node_id], fh[node_id], start_state,
                             th.branches[node_id], gen, proposed_path[node_id],
                             proposal);
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

/* Counterpart of downward_sampling_branch */
static double
proposal_prob_branch(const vector<SegmentInfo> &seg_info,
                     const FelsHelper &fh, const Path &path,
                     const size_t proposal) {

  const size_t n_intervals = seg_info.size();

  size_t start_state = path.init_state, end_state = path.init_state;
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
    const CTMarkovModel ctmm(seg_info[i].rate0, seg_info[i].rate1);
    
    /* proposal: 0 - direct, 1 - unif, 2 - poisson, 3 - forward */
    double interval_prob;
    if (proposal == 2) {
      interval_prob = // Poisson
      end_cond_sample_Poisson_prob(ctmm, path.jumps, start_state, end_state,
                                   start_time, end_time,
                                   start_jump, end_jump);
    } else if (proposal == 3) {
      interval_prob = // Forward
      forward_sample_prob(ctmm, path.jumps, start_state, end_state,
                          start_time, end_time, start_jump, end_jump);
    } else {
      interval_prob = // Direct
      end_cond_sample_prob(ctmm, path.jumps, start_state, end_state,
                           start_time, end_time, start_jump, end_jump);
    }

    log_prob += interval_prob;
    assert(std::isfinite(log_prob));

    vector<vector<double> > P;
    ctmm.get_trans_prob_mat(seg_info[i].len, P);

    // p0 = P_v(j, k) x q_k(v)/p_j(v) [along a branch, q[i]=p[i+1]
    const double p0 =
      P[start_state][0]/fh.p[i][start_state]*((i == n_intervals-1) ?
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


/* Compute proposal probabilities of uniformization algorithm */
static double
proposal_prob_branch_unif(const vector<SegmentInfo> &seg_info,
                          const FelsHelper &fh, const Path &path) {
  
  const size_t n_intervals = seg_info.size();
  
  size_t start_state = path.init_state, end_state = path.init_state;
  double start_time = 0.0, end_time = 0.0;
  size_t start_jump = 0, end_jump = 0;
  
  double log_prob = 0.0;
  
  vector<mixJump> mixtype_jumps(path.mjumps);
  if (path.jumps.size() > 0 && path.mjumps.size() == 0) {
    // This is an initial path without any virtual jump information
    for (size_t jump_id = 0; jump_id < path.jumps.size(); jump_id++)
      mixtype_jumps.push_back(mixJump(true, path.jumps[jump_id]));
  }

  for (size_t i = 0; i < n_intervals; ++i) {
    size_t n_real_jumps = 0;

    // get the end_time by adding the segment length to the start_time
    end_time += seg_info[i].len;
    
    // get the end_jump as the first jump occurring after end_time
    while (end_jump < mixtype_jumps.size() &&
           mixtype_jumps[end_jump].time < end_time) {
      n_real_jumps += (mixtype_jumps[end_jump].type);
      ++end_jump;
    }
    
    if (n_real_jumps % 2 == 1)
      end_state = complement_state(end_state);
    
    // calculate the probability for the end-conditioned path
    const CTMarkovModel ctmm(seg_info[i].rate0, seg_info[i].rate1);
    const double interval_prob =
    log(end_cond_sample_unif_prob(ctmm, mixtype_jumps, start_state, end_state,
                              start_time, end_time, start_jump, end_jump));
    log_prob += interval_prob;
    assert(std::isfinite(log_prob));
    
    vector<vector<double> > P;
    ctmm.get_trans_prob_mat(seg_info[i].len, P);
    
    // p0 = P_v(j, k) x q_k(v)/p_j(v) [along a branch, q[i]=p[i+1]
    const double p0 =
    P[start_state][0]/fh.p[i][start_state]*((i == n_intervals-1) ?
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
              const TreeHelper &th,
              const size_t site_id,
              const vector<vector<Path> > &paths,
              const vector<vector<double> > &horiz_trans_prob,
              const vector<FelsHelper> &fh,
              const vector<vector<SegmentInfo> > &seg_info,
              const vector<Path> &the_path, const size_t proposal) {

  // compute posterior probability of state 0 at root node
  const double root_p0 =
    root_post_prob0(site_id, paths[1], horiz_trans_prob, fh[0].q);
  
  double log_prob = the_path[1].init_state ? log(1.0 - root_p0) : log(root_p0);
  
  // process the paths above each node (except the root)
  for (size_t node_id = 1; node_id < th.n_nodes; ++node_id) {
    /* proposal: 0 - direct, 1 - unif, 2 - poisson, 3 - forward */
    if (proposal == 1)
      log_prob += proposal_prob_branch_unif(seg_info[node_id], fh[node_id],
                                        the_path[node_id]);
    else
      log_prob += proposal_prob_branch(seg_info[node_id], fh[node_id],
                                   the_path[node_id], proposal);
  }
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
log_accept_rate(const EpiEvoModel &mod, const TreeHelper &th,
                const size_t site_id,
                const vector<vector<Path> > &paths,
                const vector<FelsHelper> &fh,
                const vector<vector<SegmentInfo> > &seg_info,
                const vector<Path> &proposed_path, const size_t proposal) {

  static const size_t n_triples = 8;

  // ADS: this is unfortunate; we need to slice/transpose the original
  // paths because of how they are organized
  vector<Path> original(th.n_nodes);
  for (size_t i = 1; i < th.n_nodes; ++i)
    original[i] = paths[i][site_id];
  const double orig_proposal =
    proposal_prob(mod.triplet_rates, th, site_id, paths,
                  mod.init_T, fh, seg_info, original, proposal);

  const double update_proposal =
    proposal_prob(mod.triplet_rates, th, site_id, paths,
                  mod.init_T, fh, seg_info, proposed_path, proposal);

  assert(site_id > 0 && site_id < paths[1].size());

  double llr = orig_proposal - update_proposal;

  if (proposal == 1 && !original[1].mjump_touched)
    llr = 0;

  /* calculate likelihood involving root states */
  const size_t rt_l = paths[1][site_id-1].init_state;
  const size_t rt_r = paths[1][site_id+1].init_state;
  const size_t rt_orig = paths[1][site_id].init_state;
  const size_t rt_prop = proposed_path[1].init_state;
  llr += (log(root_prior_lh(rt_l, rt_prop, rt_r, mod.init_T)) -
  log(root_prior_lh(rt_l, rt_orig, rt_r, mod.init_T)));

  /* calculate likelihood involving internal intervals */
  vector<double> D_orig(n_triples), J_orig(n_triples);
  vector<double> D_prop(n_triples), J_prop(n_triples);
  vector<double> scaled_rates(n_triples);
  for (size_t i = 1; i < th.n_nodes; ++i) {
    // sufficient stats for current pentet using original path at mid
    vector<Path>::const_iterator opth(paths[i].begin() + site_id);
    fill_n(begin(D_orig), n_triples, 0.0);
    fill_n(begin(J_orig), n_triples, 0.0);
    if (site_id > 1)
      add_sufficient_statistics(*(opth - 2), *(opth - 1), *opth, J_orig, D_orig);
    add_sufficient_statistics(*(opth - 1), *opth, *(opth + 1), J_orig, D_orig);
    if (site_id < paths[i].size() - 2)
      add_sufficient_statistics(*opth, *(opth + 1), *(opth + 2), J_orig, D_orig);

    // sufficient stats for current pentet using proposed path at mid
    vector<Path>::const_iterator ppth(proposed_path.begin() + i);
    fill_n(begin(D_prop), n_triples, 0.0);
    fill_n(begin(J_prop), n_triples, 0.0);
    if (site_id > 1)
      add_sufficient_statistics(*(opth - 2), *(opth - 1), *ppth, J_prop, D_prop);
    add_sufficient_statistics(*(opth - 1), *ppth, *(opth + 1), J_prop, D_prop);
    if (site_id < paths[i].size() - 2)
      add_sufficient_statistics(*ppth, *(opth + 1), *(opth + 2), J_prop, D_prop);
    
    // add difference in log-likelihood for the proposed vs. original
    // to the Hastings ratio
    llr += (log_likelihood(mod.triplet_rates, J_prop, D_prop) -
            log_likelihood(mod.triplet_rates, J_orig, D_orig));
  }
  return llr;
}

////////////////////////////////////////////////////////////////////////////////
///////////////          single MCMC iteration         /////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* Metropolis-Hastings sampling at single site */
bool
Metropolis_Hastings_site(const EpiEvoModel &the_model, const TreeHelper &th,
                         const size_t site_id, vector<vector<Path> > &paths,
                         std::mt19937 &gen, vector<Path> &proposed_path,
                         const size_t proposal, const bool FIX_ROOT) {
  // get rates and lengths each interval [seg_info: node x site]
  vector<vector<SegmentInfo> > seg_info(th.n_nodes);
  for (size_t node_id = 1; node_id < th.n_nodes; ++node_id) {
    collect_segment_info(the_model.triplet_rates,
                         paths[node_id][site_id - 1],
                         paths[node_id][site_id + 1], seg_info[node_id]);
  }


  // upward pruning and downward sampling [fh: one for each node]
  vector<FelsHelper> fh;
  pruning(th, site_id, paths, seg_info, fh);
  downward_sampling(th, site_id, paths, the_model.init_T, seg_info, fh, gen,
                    proposed_path, proposal, FIX_ROOT);

  // acceptance rate
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  const double u = unif(gen) ;
  const double log_acc_rate =
  log_accept_rate(the_model, th, site_id, paths, fh, seg_info, proposed_path,
                  proposal);
  
  bool accepted = false;
  if (log_acc_rate >= 0 || u < exp(log_acc_rate))
    accepted = true;

  // if accepted, replace old path with proposed one.
  if (accepted)
    for (size_t i = 1; i < th.n_nodes; ++i)
      paths[i][site_id] = proposed_path[i];

  return accepted;
}
