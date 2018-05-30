/* Copyright (C) 2018 University of Southern California
 *                    Jianghan Qu and Andrew D Smith
 *
 * Author: Andrew D. Smith and Jianghan Qu
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
#include <algorithm>   // std::lower_bound,
#include <iostream>
#include <cmath>
#include <random>
#include <functional>

#include "SingleSampler.hpp"
#include "PhyloTreePreorder.hpp"
#include "Path.hpp"

#include "StateSeq.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::pair;
using std::make_pair;

// collect rates and interval lengths
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
process_branch_above(const vector<SegmentInfo> seg_info, FelsHelper &fh) {

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
pruning_branch(const vector<double> &triplet_rates, const TreeHelper &th,
               const size_t site_id, const size_t node_id,
               const vector<Path> &paths,
               const vector<SegmentInfo> &seg_info,
               vector<FelsHelper> &fh) {

  assert(fh.size() == th.n_nodes);

  /* first calculate q */
  vector<double> q(2, 1.0);
  if (is_leaf(th.subtree_sizes[node_id])) {
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
  if (!is_root(node_id))
    process_branch_above(seg_info, fh[node_id]);
}

/* Iterates over all nodes, in reverse pre-order, and for each it sets
 * the "p" and "q" from Felsenstein's algorithm, requiring for each
 * node that the children have already been processed (hence reversing
 * the order).
 */
void
pruning(const vector<double> &triplet_rates,
        const TreeHelper &th,
        const size_t site_id,
        const vector<vector<Path> > &paths,
        const vector<vector<SegmentInfo> > &seg_info,
        vector<FelsHelper> &fh) {

  // avoid recursion by iterating backwards ensuring all child nodes
  // are processed before their parents
  fh.resize(th.n_nodes);
  for (size_t i = th.n_nodes; i > 0; --i) {
    const size_t node_id = i - 1;
    pruning_branch(triplet_rates, th, site_id, node_id,
                   paths[node_id], seg_info[node_id], fh);
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


// downward sampling on a single branch, given the initial state
static void
downward_sampling_branch(const vector<SegmentInfo> &seg_info,
                         const FelsHelper &fh,
                         const size_t start_state,
                         const double branch_length,
                         std::mt19937 &gen,
                         Path &sampled_path) {

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

    const double p0 = P[prev_state][0]/fh.p[i][prev_state]*
      ((i == n_intervals - 1) ? fh.q[0] : fh.p[i + 1][0]);

    const size_t sampled_state = (unif(gen) > p0);

    const CTMarkovModel ctmm(seg_info[i].rate0, seg_info[i].rate1);
    end_cond_sample(ctmm, prev_state, sampled_state, seg_info[i].len, gen,
                    sampled_path.jumps, time_passed);

    // prepare for next interval
    time_passed += seg_info[i].len;
    prev_state = sampled_state;
  }
}

void
downward_sampling(const vector<double> &triplet_rates,
                  const TreeHelper &th,
                  const size_t site_id,
                  const vector<vector<Path> > &paths,
                  const vector<vector<double> > &horiz_trans_prob,
                  const vector<vector<SegmentInfo> > &seg_info,
                  const vector<FelsHelper> &fh,
                  std::mt19937 &gen,
                  vector<Path> &proposed_path) {

  // compute posterior probability at root node
  const double root_p0 =
    root_post_prob0(site_id, paths[0], horiz_trans_prob, fh[0].q);

  proposed_path = vector<Path>(th.n_nodes);
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  proposed_path.front() = Path(unif(gen) > root_p0, 0.0); // root

  for (size_t node_id = 1; node_id < th.n_nodes; ++node_id) {
    const size_t start_state = proposed_path[th.parent_ids[node_id]].end_state();
    downward_sampling_branch(seg_info[node_id], fh[node_id], start_state,
                             th.branches[node_id], gen, proposed_path[node_id]);
  }
}


////////////////////////////////////////////////////////////////////////////////
///////////////        Compute acceptance rate         /////////////////////////
////////////////////////////////////////////////////////////////////////////////

/*counterpart of downward_sampling_branch*/
static double
proposal_prob_branch(const vector<SegmentInfo> &seg_info,
                     const FelsHelper &fh, const Path &path) {

  const size_t n_intervals = seg_info.size();

  size_t a = path.init_state; // state at top of branch
  double time_passed = 0;
  size_t start_jump = 0;

  double prob = 1.0;
  for (size_t i = 0; i < n_intervals; ++i) {

    // find the range of jumps to evaluate using the current
    // interval's rate and duration
    size_t end_jump = start_jump;
    while (end_jump < path.jumps.size() - 1 &&
           path.jumps[end_jump + 1] < time_passed + seg_info[i].len)
      ++end_jump;

    const size_t b = (end_jump - start_jump) % 2 == 0 ? a : complement_state(a);

    const CTMarkovModel ctmm(seg_info[i].rate0, seg_info[i].rate1);
    prob *= end_cond_sample_prob(ctmm, a, b, time_passed + seg_info[i].len,
                                 path.jumps, start_jump, end_jump, time_passed);

    vector<vector<double> > P;
    ctmm.get_trans_prob_mat(seg_info[i].len, P);

    // p0 = P_v(j, k) x q_k(v)/p_j(v) [along a branch, q[i]=p[i+1]
    const double p0 = P[a][0]/fh.p[i][a]*((i == n_intervals-1) ? fh.q[0] : fh.p[i+1][0]);

    prob *= (b == 0) ? p0 : 1.0 - p0;

    // prepare for next interval
    start_jump = end_jump + 1;
    time_passed += seg_info[i].len;
    a = b;
  }
  return prob;
}


// compute proposal prob
static double
proposal_prob(const vector<double> &triplet_rates,
              const TreeHelper &th,
              const size_t site_id,
              const vector<vector<Path> > &paths,
              const vector<vector<double> > &horiz_trans_prob,
              const vector<FelsHelper> &fh,
              const vector<vector<SegmentInfo> > &seg_info,
              const vector<Path> &the_path) {

  // compute posterior probability of state 0 at root node
  const double root_p0 =
    root_post_prob0(site_id, paths[0], horiz_trans_prob, fh[0].q);

  double prob = the_path[0].init_state ? 1.0 - root_p0 : root_p0;

  // process the paths above each node (except the root)
  for (size_t node_id = 1; node_id < th.n_nodes; ++node_id)
    prob *= proposal_prob_branch(seg_info[node_id], fh[node_id], the_path[node_id]);
  return prob;
}


static double
log_lik_ratio(const vector<double> &rates,
              const PathContextStat &pcs_num,
              const PathContextStat &pcs_denom) {
  static const size_t n_triples = 8;
  double result = 0.0;
  for (size_t i = 0; i < n_triples; ++i) {
    const double J_diff = pcs_num.jumps_in_context[i] - pcs_denom.jumps_in_context[i];
    const double D_diff = pcs_num.time_in_context[i] - pcs_denom.time_in_context[i];
    result += J_diff*log(rates[i]) - D_diff*rates[i];
  }
  return result;
}


// compute acceptance rate
double
log_accept_rate(const EpiEvoModel &the_model, const TreeHelper &th,
                const size_t site_id,
                const vector<vector<Path> > &paths,
                const vector<FelsHelper> &fh,
                const vector<vector<SegmentInfo> > &seg_info,
                const vector<Path> &proposed_path) {

  vector<Path> original(th.n_nodes);
  for (size_t i = 0; i < th.n_nodes; ++i)
    original[i] = paths[i][site_id];

  const double prob_old =
    proposal_prob(the_model.triplet_rates, th,
                  site_id, paths, the_model.init_T, fh, seg_info, original);

  const double prob_new =
    proposal_prob(the_model.triplet_rates, th,
                  site_id, paths, the_model.init_T, fh, seg_info, proposed_path);

  double lr = log(prob_old) - log(prob_new);

  for (size_t i = 1; i < th.n_nodes; ++i) {
    const Path l =  paths[i][site_id - 1];
    const Path ll = paths[i][site_id - 2];
    const Path r =  paths[i][site_id + 1];
    const Path rr = paths[i][site_id + 2];
    PathContextStat pcs_old(l, paths[i][site_id], r);
    PathContextStat pcs_new(l, proposed_path[i], r);
    PathContextStat pcs_old_l(ll, l, paths[i][site_id]);
    PathContextStat pcs_new_l(ll, l, proposed_path[i]);
    PathContextStat pcs_old_r(paths[i][site_id], r, rr);
    PathContextStat pcs_new_r(proposed_path[i], r, rr);
    lr +=
      log_lik_ratio(the_model.triplet_rates, pcs_new, pcs_old) +
      log_lik_ratio(the_model.triplet_rates, pcs_new_l, pcs_old_l) +
      log_lik_ratio(the_model.triplet_rates, pcs_new_r, pcs_old_r);
  }

  return lr;
}


void
gibbs_site(const EpiEvoModel &the_model, const TreeHelper &th,
           const size_t site_id, vector<vector<Path> > &paths,
           std::mt19937 &gen, vector<Path> &proposed_path) {

  // get rates and lengths each interval [seg_info: node x site]
  vector<vector<SegmentInfo> > seg_info(th.n_nodes);
  for (size_t node_id = 1; node_id < th.n_nodes; ++node_id)
    collect_segment_info(the_model.triplet_rates,
                         paths[node_id][site_id - 1],
                         paths[node_id][site_id + 1], seg_info[node_id]);

  // upward pruning and downward sampling [fh: one for each node]
  vector<FelsHelper> fh;
  pruning(the_model.triplet_rates, th, site_id, paths, seg_info, fh);

  downward_sampling(the_model.triplet_rates, th, site_id,
                    paths, the_model.init_T, seg_info, fh, gen, proposed_path);

  // acceptance rate
  const double log_acc_rate =
    log_accept_rate(the_model, th, site_id, paths, fh, seg_info, proposed_path);

  std::uniform_real_distribution<double> unif(0.0, 1.0);

  if (log_acc_rate >= 0 || unif(gen) < exp(log_acc_rate))
    for (size_t i = 0; i < th.n_nodes; ++i)
      paths[i][site_id] = proposed_path[i];
}
