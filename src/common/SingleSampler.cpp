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


////////////////////////////////////////////////////////////////////////////////
///////////////////           UPWARD PRUNING           /////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* rates are 2 muation rates, for patterns L0R and L1R. The "q" and
 * "p" values are as defined for Felsenstein's algorithm generally.
 */
static void
upward(const pair<double, double> &rates, const double time_interval,
       const vector<double> &q, vector<double> &p) {

  assert(time_interval > 0);

  vector<vector<double> > P; // transition matrix
  continuous_time_trans_prob_mat(rates.first, rates.second, time_interval, P);

  // p <- P*q
  p.resize(2);
  p[0] = P[0][0]*q[0] + P[0][1]*q[1];
  p[1] = P[1][0]*q[0] + P[1][1]*q[1];
}


// collect rates and interval lengths
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

// all_p: branch x interval x state
// for a single branch
void
pruning_branch(const vector<double> &triplet_rates,
               const vector<size_t> &subtree_sizes,
               const size_t site_id,
               const size_t node_id,
               const vector<vector<Path> > &all_paths,
               const vector<interval_rate_pairs> &all_interval_rates,
               const vector<vector<double> > &all_interval_lengths,
               vector<vector<vector<double> > > &p) {

  assert(p.size() == subtree_sizes.size());

  vector<double> q(2, 1.0);
  if (is_leaf(subtree_sizes[node_id])) {
    const bool leaf_state = all_paths[node_id][site_id].end_state();
    q[0] = (leaf_state == false) ? 1.0 : 0.0;
    q[1] = (leaf_state == true)  ? 1.0 : 0.0;
  }
  else {

    for (size_t ch_id = node_id + 1; ch_id < subtree_sizes[node_id];
         ch_id += subtree_sizes[ch_id]) {
      pruning_branch(triplet_rates, subtree_sizes, site_id, ch_id,
                     all_paths, all_interval_rates, all_interval_lengths, p);
      q[0] *= p[ch_id].back()[0];
      q[1] *= p[ch_id].back()[1];
    }
  }

  if (!is_root(node_id)) {
    const size_t n_intervals = all_interval_lengths[node_id].size();
    assert(node_id == 0 || n_intervals > 0);
    p[node_id] = vector<vector<double> >(n_intervals, {0.0, 0.0});

    // directly use q for final interval (back)
    upward(all_interval_rates[node_id].back(),
           all_interval_lengths[node_id].back(), q, p[node_id].back());

    // iterate backwards (upwards) with: q[i-1] = p[i]
    for (size_t i = n_intervals - 1; i > 0; --i)
      upward(all_interval_rates[node_id][i-1],
             all_interval_lengths[node_id][i-1], p[node_id][i], p[node_id][i-1]);
  }
}


////////////////////////////////////////////////////////////////////////////////
///////////////            DOWNWARD SAMPLING        ////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// downward sampling on a single branch, given the initial state
void
downward_sampling_branch(const interval_rate_pairs &interval_rates,
                         const vector<double> &interval_lengths,
                         const vector<size_t> &subtree_sizes,
                         const size_t site,
                         const size_t node_id,
                         const vector<vector<Path> > &all_paths,
                         const vector<vector<vector<double> > > &all_p,
                         std::mt19937 &gen,
                         vector<Path> &new_path) {

  const size_t n_intervals = interval_rates.size();

  size_t parent_state = new_path[node_id].init_state;
  double time_passed = 0;

  // end-conditioned sampling helpers
  for (size_t m = 0; m < n_intervals - 1; ++m) {

    // compute conditional posterior probability
    vector<vector<double> > P; // transition prob matrix
    continuous_time_trans_prob_mat(interval_rates[m].first,
                                   interval_rates[m].second,
                                   interval_lengths[m], P);
    double p0 = (all_p[node_id][m+1][0] * P[parent_state][0] /
                 all_p[node_id][m][parent_state]);

    // generate random state at break point
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    bool new_state = (unif(gen) > p0);

    // prepare helper values
    const CTMarkovModel ctmm(interval_rates[m]);

    // generate path
    vector<double> jump_times;
    end_cond_sample(ctmm, parent_state, (size_t)(new_state),
                    interval_lengths[m], gen, jump_times);

    // append jump_times to new_path
    for (size_t i = 0; i < jump_times.size(); ++i)
      new_path[node_id].jumps.push_back(time_passed + jump_times[i]);

    // prepare for next interval
    time_passed += interval_lengths[m];
    parent_state = new_state;
  }

  if (is_leaf(subtree_sizes[node_id])) {
    const size_t leaf_state = all_paths[node_id][site].end_state();
    // generate path
    vector<double> jump_times;
    const size_t m = n_intervals - 1;
    // prepare helper values
    const CTMarkovModel ctmm(interval_rates[m]);
    end_cond_sample(ctmm, parent_state, leaf_state,
                    interval_lengths[m], gen, jump_times);
    // append jump_times to new_path
    for (size_t i = 0; i < jump_times.size(); ++i)
      new_path[node_id].jumps.push_back(time_passed + jump_times[i]);

    assert(new_path[node_id].end_state() == leaf_state);
  }
  else {
    // last interval requires information from children
    vector<size_t> children;
    get_children(node_id, subtree_sizes, children);
    const size_t m = n_intervals - 1;
    // compute conditional posterior probability
    vector<vector<double> > P; // transition matrix
    continuous_time_trans_prob_mat(interval_rates[m].first,
                                   interval_rates[m].second,
                                   interval_lengths[m], P);
    double p0 = 1.0;
    for (size_t idx = 0; idx < children.size(); ++idx)
      p0 *= all_p[children[idx]][0][0];
    p0 *= P[parent_state][0] / all_p[node_id][m][parent_state];

    assert(p0 > 0 && p0 < 1);

    // generate random state at break point
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    bool new_state = (unif(gen) > p0);
    // set new initial states for children branches
    for (size_t idx = 0; idx < children.size(); ++idx) {
      new_path[children[idx]].init_state = new_state;
      new_path[children[idx]].tot_time =
        all_paths[children[idx]][site - 1].tot_time;
    }

    // prepare helper values
    const CTMarkovModel ctmm(interval_rates[m]);
    // generate path
    vector<double> jump_times;
    end_cond_sample(ctmm, parent_state, (size_t)(new_state),
                    interval_lengths[m], gen, jump_times);
    // append jump_times to new_path
    for (size_t i = 0; i < jump_times.size(); ++i) {
      new_path[node_id].jumps.push_back(time_passed + jump_times[i]) ;
    }
  }
}


/* Compute posterior probability of state 0 at root node (i.e. the
 * init_state of any/all children) for particular site using
 * information from the upward stage of Felsenstein's algorithm.
 */
double
root_post_prob0(const vector<size_t> &children,
                const size_t site_id,
                const vector<Path> &the_paths,
                const vector<vector<double> > &horiz_trans_prob,
                const vector<vector<vector<double> > > &all_p) {

  const size_t left_state = the_paths[site_id - 1].init_state;
  const size_t right_state = the_paths[site_id + 1].init_state;

  double p0 = (horiz_trans_prob[left_state][0]*
               horiz_trans_prob[0][right_state]);
  double p1 = (horiz_trans_prob[left_state][1]*
               horiz_trans_prob[1][right_state]);

  for (size_t idx = 0; idx < children.size(); ++idx) {
    p0 *= all_p[children[idx]][0][0];
    p1 *= all_p[children[idx]][0][1];
  }

  return p0/(p0 + p1);
}

void
downward_sampling(const vector<double> &triplet_rates,
                  const vector<size_t> &subtree_sizes,
                  const size_t site_id,
                  const vector<vector<Path> > &all_paths,
                  const vector<vector<double> > &horiz_trans_prob,
                  const vector<vector<vector<double> > > &all_p,
                  std::mt19937 &gen,
                  vector<Path> &proposed_path) {

  proposed_path = vector<Path>(subtree_sizes.size(), Path());

  const size_t root_id = 0; //all_paths[0] is empty
  vector<size_t> children;
  get_children(root_id, subtree_sizes, children);

  // compute posterior probability at root node
  const double root_p0 = root_post_prob0(children, site_id, all_paths[1],
                                         horiz_trans_prob, all_p);

  // sample new root state
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  const bool new_root_state = (unif(gen) > root_p0);

  for (size_t idx = 0; idx < children.size(); ++idx) {
    proposed_path[children[idx]].init_state = new_root_state;
    proposed_path[children[idx]].tot_time = all_paths[children[idx]][site_id - 1].tot_time;
  }

  // preorder traversal of the tree
  for (size_t node_id = 1; node_id < subtree_sizes.size(); ++node_id) {
    // collect relevant transition rates for each interval
    interval_rate_pairs interval_rates; //interval x transition type
    vector<double> interval_lengths;
    collect_interval_rates_and_lengths(triplet_rates, all_paths[node_id][site_id-1],
                    all_paths[node_id][site_id+1],
                    interval_rates, interval_lengths);
    downward_sampling_branch(interval_rates, interval_lengths, subtree_sizes,
                             site_id, node_id, all_paths, all_p,
                             gen, proposed_path);
  }
}


////////////////////////////////////////////////////////////////////////////////
///////////////        Compute acceptance rate         /////////////////////////
////////////////////////////////////////////////////////////////////////////////


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
                                 jumps, start_jump, end_jump, time_passed);

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
              const vector<size_t> &subtree_sizes,
              const size_t site_id,
              const vector<vector<Path> > &all_paths,
              const vector<vector<double> > &horiz_trans_prob,
              const vector<vector<vector<double> > > &all_p,
              const vector<Path> &proposed_path,
              double &prob_old, double &prob_new) {

  static const size_t root_id = 0;

  const size_t n_nodes = subtree_sizes.size();

  vector<size_t> children;
  get_children(root_id, subtree_sizes, children);

  // compute posterior probability of state 0 at root node
  const double root_p0 = root_post_prob0(children, site_id, all_paths[1],
                                         horiz_trans_prob, all_p);

  prob_old = all_paths[1][site_id].init_state ? 1.0 - root_p0 : root_p0;
  prob_new = proposed_path[1].init_state      ? 1.0 - root_p0 : root_p0;

  // process the paths above each node (except the root)
  for (size_t node_id = 1; node_id < n_nodes; ++node_id) {

    // collect transition rates and lengths for each interval
    interval_rate_pairs interval_rates;
    vector<double> interval_lengths;
    collect_interval_rates_and_lengths(triplet_rates,
                                       all_paths[node_id][site_id - 1],
                                       all_paths[node_id][site_id + 1],
                                       interval_rates, interval_lengths);

    proposal_prob_branch(interval_rates, interval_lengths,
                         subtree_sizes, site_id, node_id, all_paths,
                         all_p, proposed_path[node_id], prob_new);

    proposal_prob_branch(interval_rates, interval_lengths,
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
                const vector<Path> &proposed_path) {

  double prob_old, prob_new;
  proposal_prob(the_model.triplet_rates, th.subtree_sizes,
                site_id, all_paths, the_model.init_T, all_p, proposed_path,
                prob_old, prob_new);

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
  vector<interval_rate_pairs> all_interval_rates(n_nodes);
  vector<vector<double> > all_interval_lengths(n_nodes);
  for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
    collect_interval_rates_and_lengths(the_model.triplet_rates,
                                       all_paths[node_id][site_id - 1],
                                       all_paths[node_id][site_id + 1],
                                       all_interval_rates[node_id],
                                       all_interval_lengths[node_id]);
  }

  // upward pruning and downward sampling
  vector<vector<vector<double> > > all_p;
  pruning_branch(the_model.triplet_rates, th.subtree_sizes, site_id, 0,
                 all_paths, all_interval_rates, all_interval_lengths, all_p);

  downward_sampling(the_model.triplet_rates, th.subtree_sizes, site_id,
                    all_paths, the_model.init_T, all_p, gen, proposed_path);

  // acceptance rate
  const double log_acc_rate =
    log_accept_rate(the_model, th, site_id, all_paths, all_p, proposed_path);

  std::uniform_real_distribution<double> unif(0.0, 1.0);

  if (log_acc_rate >= 0 || unif(gen) < exp(log_acc_rate))
    for (size_t i = 1; i < th.subtree_sizes.size(); ++i)
      all_paths[i][site_id] = proposed_path[i];
}
