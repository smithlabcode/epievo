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
// Liz: where is interval_rate_pairs defined?
void
collect_interval_rates_and_lengths(const vector<double> &triplet_rates,
                                   const Path &l, const Path &r,
                                   interval_rate_pairs &interval_rates,
                                   vector<double> &interval_lengths) {
  Environment env(l, r);
  const size_t n_intervals = env.left.size();
  interval_rates = interval_rate_pairs(n_intervals, make_pair(0.0, 0.0));
  interval_lengths = vector<double>(n_intervals, 0.0);

  for (size_t i = 0; i < n_intervals; ++i) {
    const size_t pattern0 = triple2idx(env.left[i], false, env.right[i]);
    const size_t pattern1 = triple2idx(env.left[i], true, env.right[i]);
    interval_rates[i].first = triplet_rates[pattern0];
    interval_rates[i].second = triplet_rates[pattern1];
    interval_lengths[i] = (i == 0) ? env.breaks[0] : env.breaks[i] - env.breaks[i-1];
    assert(interval_lengths[i] > 0);
  }
}


////////////////////////////////////////////////////////////////////////////////
///////////////////           UPWARD PRUNING           /////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* rates are 2 muation rates, for patterns L0R and L1R. The "q" and
 * "p" values are as defined for Felsenstein's algorithm generally.
 */
static void
process_interval(const pair<double, double> &rates, const double time_interval,
                 const vector<double> &q, vector<double> &p) {

  assert(time_interval > 0);

  vector<vector<double> > P; // transition matrix
  continuous_time_trans_prob_mat(rates.first, rates.second, time_interval, P);

  // p <- P*q
  p.resize(2);
  p[0] = P[0][0]*q[0] + P[0][1]*q[1];
  p[1] = P[1][0]*q[0] + P[1][1]*q[1];
}


static void
process_branch_above(const interval_rate_pairs interval_rates,
                     const vector<double> interval_lengths,
                     const vector<double> &q, vector<vector<double> > &p) {

  const size_t n_intervals = interval_lengths.size();
  p = vector<vector<double> >(n_intervals, {0.0, 0.0});
  assert(n_intervals > 0);

  // directly use q for final interval (back of p[node_id])
  process_interval(interval_rates.back(), interval_lengths.back(), q, p.back());

  // iterate backwards/upwards with: q[i-1] = p[i]
  for (size_t i = n_intervals - 1; i > 0; --i)
    process_interval(interval_rates[i-1], interval_lengths[i-1], p[i], p[i-1]);
}


/* p: node x interval x state
 * q: state
 */
void
pruning_branch(const vector<double> &triplet_rates,
               const TreeHelper &th,
               const size_t site_id,
               const size_t node_id,
               const vector<Path> &paths,
               const interval_rate_pairs &interval_rates,
               const vector<double> &interval_lengths,
               vector<vector<vector<double> > > &p,
               vector<double> &q) {

  assert(p.size() == th.n_nodes && q.size() == th.n_nodes);

  /* first calculate q */
  q = vector<double>(2, 1.0);
  if (is_leaf(th.subtree_sizes[node_id])) {
    const bool leaf_state = paths[site_id].end_state();
    q[0] = (leaf_state == false) ? 1.0 : 0.0;
    q[1] = (leaf_state == true)  ? 1.0 : 0.0;
  }
  else {
    for (ChildSet ch_id(th.subtree_sizes, node_id); ch_id.good(); ++ch_id) {
      q[0] *= p[*ch_id].back()[0];
      q[1] *= p[*ch_id].back()[1];
    }
  }

  /* now calculate p */
  if (!is_root(node_id))
    process_branch_above(interval_rates, interval_lengths, q, p[node_id]);
}


// all_p: branch x interval x state
// Liz: (to-do) keep q values for each node
void
pruning(const vector<double> &triplet_rates,
        const TreeHelper &th,
        const size_t site_id,
        const vector<vector<Path> > &all_paths,
        const vector<interval_rate_pairs> &interval_rates,
        const vector<vector<double> > & interval_lengths,
        vector<vector<vector<double> > > &p,
        vector<vector<double> > &q) {

  p.resize(th.n_nodes);
  q.resize(th.n_nodes);

  // avoid recursion by iterating backwards ensuring all child nodes
  // are processed before their parents
  for (size_t i = th.n_nodes; i > 0; --i) {
    const size_t node_id = i - 1;
    pruning_branch(triplet_rates, th, site_id, node_id,
                   all_paths[node_id], interval_rates[node_id],
                   interval_lengths[node_id], p, q[node_id]);
  }
}


////////////////////////////////////////////////////////////////////////////////
///////////////            DOWNWARD SAMPLING        ////////////////////////////
////////////////////////////////////////////////////////////////////////////////


/* Compute posterior probability of state 0 at root node (i.e. the
 * init_state of any/all children) for particular site using
 * information from the upward stage of Felsenstein's algorithm.
 */
double
root_post_prob0(const size_t site_id,
                const vector<Path> &the_paths,
                const vector<vector<double> > &horiz_trans_prob,
                const vector<double> &q) {

  const size_t left_state = the_paths[site_id - 1].init_state;
  const size_t right_state = the_paths[site_id + 1].init_state;

  const double p0 =
    (horiz_trans_prob[left_state][0]*horiz_trans_prob[0][right_state])*q[0];
  const double p1 =
    (horiz_trans_prob[left_state][1]*horiz_trans_prob[1][right_state])*q[1];
  return p0/(p0 + p1);
}


// downward sampling on a single branch, given the initial state
static void
downward_sampling_branch(const interval_rate_pairs &interval_rates,
                         const vector<double> &interval_lengths,
                         const vector<vector<double> > &p,
                         const vector<double> &q,
                         const size_t start_state,
                         const double branch_length,
                         std::mt19937 &gen,
                         Path &sampled_path) {

  const size_t n_intervals = interval_rates.size();

  sampled_path.init_state = start_state;
  sampled_path.tot_time = branch_length;

  std::uniform_real_distribution<double> unif(0.0, 1.0);

  size_t prev_state = sampled_path.init_state;
  double time_passed = 0;
  for (size_t m = 0; m < n_intervals; ++m) {

    vector<vector<double> > P;
    continuous_time_trans_prob_mat(interval_rates[m].first,
                                   interval_rates[m].second, interval_lengths[m], P);

    const double p0 = P[prev_state][0]/p[m][prev_state]*
      ((m == n_intervals - 1) ? q[0] : p[m + 1][0]);
    const size_t sampled_state = (unif(gen) > p0);

    const CTMarkovModel ctmm(interval_rates[m]);
    end_cond_sample(ctmm, prev_state, sampled_state, interval_lengths[m], gen,
                    sampled_path.jumps, time_passed);

    // prepare for next interval
    time_passed += interval_lengths[m];
    prev_state = sampled_state;
  }
}

void
downward_sampling(const vector<double> &triplet_rates,
                  const TreeHelper &th,
                  const size_t site_id,
                  const vector<vector<Path> > &all_paths,
                  const vector<vector<double> > &horiz_trans_prob,
                  const vector<interval_rate_pairs> &interval_rates,
                  const vector<vector<double> > &interval_lengths,
                  const vector<vector<vector<double> > > &p,
                  const vector<vector<double> > &q,
                  std::mt19937 &gen,
                  vector<Path> &proposed_path) {

  proposed_path = vector<Path>(th.n_nodes);

  // compute posterior probability at root node
  const double root_p0 =
    root_post_prob0(site_id, all_paths[0], horiz_trans_prob, q[node_id]);
  std::uniform_real_distribution<double> unif(0.0, 1.0);

  proposed_path.front() = Path(unif(gen) > root_p0, 0.0); // root

  for (size_t node_id = 1; node_id < th.n_nodes; ++node_id) {
    const size_t start_state = proposed_path[th.parent_ids[node_id]].end_state();
    downward_sampling_branch(interval_rates[node_id],
                             interval_lengths[node_id],
                             p[node_id], q[node_id],
                             start_state, th.branches[node_id],
                             gen, proposed_path);
  }
}


////////////////////////////////////////////////////////////////////////////////
///////////////        Compute acceptance rate         /////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* Liz: to delete
static void
select_jumps(const vector<double> jumps,
             const double t0, const double t1,
             vector<double> & subset_jumps) {
  // jumps is a vector of ascending positive values

  vector<double>::const_iterator low, up;
  low = std::upper_bound(jumps.begin(), jumps.end(), t0); // first > t0
  up = std::lower_bound(jumps.begin(), jumps.end(), t1); // first >= t1
  subset_jumps = vector<double>(low, up);
  // change reference point to t0
  for (size_t i = 0; i < subset_jumps.size(); ++i)
    subset_jumps[i] -= t0;
}
*/

/*counterpart of downward_sampling_branch*/
void
proposal_prob_branch(const interval_rate_pairs &interval_rates,
                     const vector<double> &interval_lengths,
                     const vector<size_t> &subtree_sizes,
                     const size_t site_id,
                     const size_t node_id,
                     const vector<vector<Path> > &all_paths,
                     const vector<vector<vector<double> > > &all_p,
                     const Path &path,
                     double &prob) {

  assert(!is_leaf(subtree_sizes[node_id]) ||
         path.end_state() == all_paths[node_id][site_id].end_state());

  const size_t n_intervals = interval_rates.size();

  size_t a = path.init_state; // state at top of branch
  double time_passed = 0;
  size_t start_jump = 0;
  for (size_t i = 0; i < n_intervals - 1; ++i) {
    // find the range of jumps to evaluate using the current
    // interval's rate and duration
    size_t end_jump = start_jump;
    while (end_jump < path.jumps.size() - 1 &&
           path.jumps[end_jump + 1] < time_passed + interval_lengths[i])
      ++end_jump;

    size_t b = (end_jump - start_jump) % 2 == 0 ? a : complement_state(a);
    const CTMarkovModel ctmm(interval_rates[i]);
    prob *= end_cond_sample_prob(ctmm, a, b, time_passed + interval_lengths[i],
                                 path.jumps, start_jump, end_jump, time_passed);

    vector<vector<double> > P;
    ctmm.get_trans_prob_mat(interval_lengths[i], P);
    // p0 = P_v(j, k) x q_k(v)/p_j(v) [along a branch, q[i]=p[i+1]
    const double p0 = P[a][0]*all_p[node_id][i + 1][0]/all_p[node_id][i][a];
    prob *= (b == 0) ? p0 : 1.0 - p0;

    // prepare for next interval
    start_jump = end_jump + 1;
    time_passed += interval_lengths[i];
    a = b;
  }

  const size_t b = path.end_state();
  const CTMarkovModel ctmm(interval_rates.back());
  prob *= end_cond_sample_prob(ctmm, a, b, interval_lengths.back(),
                               start_jump, jumps.size() - 1, jumps, time_passed);

  if (!is_leaf(subtree_sizes[node_id])) {
    vector<vector<double> > P;
    ctmm.get_trans_prob_mat(interval_lengths.back(), P);

    double p0 = 1.0;
    for (auto c = node_id + 1; c < subtree_sizes[node_id]; c += subtree_sizes[c])
      p0 *= all_p[c].front()[0];
    p0 *= P[a][0]/all_p[node_id].back()[a];

    prob *= (b == 0) ? p0 : 1.0 - p0;
  }
}


// compute proposal prob
void
proposal_prob(const vector<double> &triplet_rates,
              const TreeHelper &th,
              const size_t site_id,
              const vector<vector<Path> > &all_paths,
              const vector<vector<double> > &horiz_trans_prob,
              const vector<vector<vector<double> > > &all_p,
              const vector<vector<double> > &q,
              const vector<interval_rate_pairs> &interval_rates,
              const vector<vector<double> > &interval_lengths,
              const vector<Path> &proposed_path,
              double &prob_old, double &prob_new) {

  // compute posterior probability of state 0 at root node
  const double root_p0 =
    root_post_prob0(site_id, all_paths[0], horiz_trans_prob, q);

  prob_old = all_paths[1][site_id].init_state ? 1.0 - root_p0 : root_p0;
  prob_new = proposed_path[1].init_state      ? 1.0 - root_p0 : root_p0;

  // process the paths above each node (except the root)
  for (size_t node_id = 1; node_id < th.n_nodes; ++node_id) {

    proposal_prob_branch(interval_rates[node_id], interval_lengths[node_id],
                         subtree_sizes, site_id, node_id, all_paths,
                         all_p, proposed_path[node_id], prob_new);

    proposal_prob_branch(interval_rates[node_id], interval_lengths[node_id],
                         subtree_sizes, site_id, node_id, all_paths,
                         all_p, all_paths[node_id][site_id], prob_old);
  }
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
                const vector<vector<Path> > &all_paths,
                const vector<vector<vector<double> > > &all_p,
                const vector<interval_rate_pairs> &interval_rates,
                const vector<vector<double> > &interval_lengths,
                const vector<Path> &proposed_path) {
  double prob_old, prob_new;
  proposal_prob(the_model.triplet_rates, th.subtree_sizes,
                site_id, all_paths, the_model.init_T, all_p,
                interval_rates[node_id], interval_lengths[node_id],
                proposed_path, prob_old, prob_new);
  double lr = log(prob_old) - log(prob_new);

  for (size_t i = 1; i < th.subtree_sizes.size(); ++i) {
    const Path l =  all_paths[i][site_id - 1];
    const Path ll = all_paths[i][site_id - 2];
    const Path r =  all_paths[i][site_id + 1];
    const Path rr = all_paths[i][site_id + 2];
    PathContextStat pcs_old(l, all_paths[i][site_id], r);
    PathContextStat pcs_new(l, proposed_path[i], r);
    PathContextStat pcs_old_l(ll, l, all_paths[i][site_id]);
    PathContextStat pcs_new_l(ll, l, proposed_path[i]);
    PathContextStat pcs_old_r(all_paths[i][site_id], r, rr);
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
           const size_t site_id,
           vector<vector<Path> > &all_paths,
           std::mt19937 &gen,
           vector<Path> &proposed_path) {

  // collect relevant transition rates for each interval
  const size_t n_nodes = th.subtree_sizes.size();
  vector<interval_rate_pairs> interval_rates(n_nodes);
  vector<vector<double> > interval_lengths(n_nodes);
  for (size_t node_id = 1; node_id < n_nodes; ++node_id)
    collect_interval_rates_and_lengths(the_model.triplet_rates,
                                       all_paths[node_id][site_id - 1],
                                       all_paths[node_id][site_id + 1],
                                       interval_rates[node_id],
                                       interval_lengths[node_id]);

  // upward pruning and downward sampling
  vector<vector<vector<double> > > p;
  vector<vector<double> > q;
  pruning_branch(the_model.triplet_rates, th.subtree_sizes, site_id, 0,
                 all_paths, interval_rates, interval_lengths, p, q);

  downward_sampling(the_model.triplet_rates, th, site_id,
                    all_paths, the_model.init_T,
                    interval_rates, interval_lengths, p, q, gen, proposed_path);

  // acceptance rate
  const double log_acc_rate =
    log_accept_rate(the_model, th, site_id, all_paths, all_p,
                    interval_rates, interval_lengths, proposed_path);

  std::uniform_real_distribution<double> unif(0.0, 1.0);

  if (log_acc_rate >= 0 || unif(gen) < exp(log_acc_rate))
    for (size_t i = 0; i < th.subtree_sizes.size(); ++i)
      all_paths[i][site_id] = proposed_path[i];
}
