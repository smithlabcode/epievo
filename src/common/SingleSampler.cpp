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

#include "SingleSampler.hpp"
#include "PhyloTreePreorder.hpp"
#include "Path.hpp"

#include "StateSeq.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;


////////////////////////////////////////////////////////////////////////////////
///////////////////           UPWARD PRUNING           /////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* rates are 2 muation rates, for patterns L0R and L1R. The "q" and
 * "p" values are as defined for Felsenstein's algorithm generally.
 */
static void
upward(const vector<double> &rates, const double time_interval,
       const vector<double> &q, vector<double> &p) {

  assert(time_interval > 0);

  vector<vector<double> > P; // transition matrix
  continuous_time_trans_prob_mat(rates[0], rates[1], time_interval, P);

  // p <- P*q
  p.resize(2);
  p[0] = P[0][0]*q[0] + P[0][1]*q[1];
  p[1] = P[1][0]*q[0] + P[1][1]*q[1];
}


// collect rates and interval lengths
void
rates_on_branch(const vector<double> &triplet_rates,
                const Path &l, const Path &r,
                vector<pair<double, double> > &interval_rates,
                vector<double> &interval_lengths) {
  Environment env(l, r);
  const size_t n_intervals = env.left.size();
  interval_rates = vector<pair<double, double> >(n_intervals, make_pair(0.0, 0.0));
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
               const size_t site,
               const size_t node_id,
               const vector<vector<Path> > &all_paths,
               const vector<vector<pair<double, double> > > &all_interval_rates,
               const vector<vector<double> > &all_interval_lengths,
               vector<vector<vector<double> > > &all_p) {

  // excluding the top end, i.e. including only the lower end of each time interval
  const size_t n_intervals = all_interval_lengths[node_id].size();
  vector<vector<double> > p(n_intervals, vector<double>(2, 0.0));

  vector<double> q(2, 0.0);

  if (is_leaf(subtree_sizes[node_id])) {  // leaf

    for (size_t i = 0; i < n_intervals; ++i) { // going upward
      const size_t interval = n_intervals - 1 - i;
      if (i == 0) { // at leaf
        const bool leaf_state = all_paths[node_id][site].end_state();
        for (size_t k = 0; k < 2; ++k) {
          q[k] = ((size_t)(leaf_state) == k) ? 1.0 : 0.0;
        }
        upward(all_interval_rates[node_id][interval],
               all_interval_lengths[node_id][interval],
               q, p[interval]);
      } else {
        // at break points within a branch: q = p[interval + 1]
        upward(all_interval_rates[node_id][interval],
               all_interval_lengths[node_id][interval],
               p[interval + 1], p[interval]);
      }
    }

    all_p[node_id].swap(p);

  } else { // internal node

    vector<size_t> children;
    get_children(node_id, subtree_sizes, children);

    // fill in all_p at children nodes
    for (size_t idx = 0; idx < children.size(); ++idx) {
      pruning_branch(triplet_rates, subtree_sizes, site,
                     children[idx], all_paths,
                     all_interval_rates, all_interval_lengths,
                     all_p);
    }

    if (node_id > 0) { // non-root internal node
      // going upward along current branch
      for (size_t i = 0; i < n_intervals; ++i) {
        const size_t interval = n_intervals - 1 - i;
        if (i == 0) {
          // at the end of the last interval, i.e. an internal species
          q = vector<double>(2, 1.0);
          for (size_t k = 0; k < 2; ++k) {
            for (size_t idx = 0; idx < children.size(); ++idx)
              q[k] *= all_p[children[idx]][0][k];
          }
          upward(all_interval_rates[node_id][interval],
                 all_interval_lengths[node_id][interval],
                 q, p[interval]);
        } else {
          // at break points within a branch: q = p[interval + 1]
          upward(all_interval_rates[node_id][interval],
                 all_interval_lengths[node_id][interval],
                 p[interval + 1], p[interval]);
        }
      }
    }
    all_p[node_id].swap(p);
  }
}



// all_p: branch x interval x state
void
pruning(const vector<double> &triplet_rates,
        const vector<size_t> &subtree_sizes,
        const size_t site,
        const vector<vector<Path> > &all_paths,
        const vector<vector<vector<double> > > &all_interval_rates,
        const vector<vector<double> > & all_interval_lengths,
        vector<vector<vector<double> > > &all_p) {

  all_p.resize(subtree_sizes.size());

  size_t child = 1;
  while(child < subtree_sizes.size()) {
    pruning_branch(triplet_rates, subtree_sizes, site, child, all_paths,
                   all_interval_rates, all_interval_lengths, all_p);
    child += subtree_sizes[child];
  }
}

////////////////////////////////////////////////////////////////////////////////
///////////////            DOWNWARD SAMPLING        ////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// downward sampling on a single branch, given the initial state
void
downward_sampling_branch(const vector<vector<double> > &interval_rates,
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
    continuous_time_trans_prob_mat(interval_rates[m][0], interval_rates[m][1],
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
    continuous_time_trans_prob_mat(interval_rates[m][0], interval_rates[m][1],
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

double
root_post_prob0(const vector<size_t> &children,
                const size_t site,
                const vector<vector<Path> > &all_paths,
                const vector<vector<double> > &root_trans_prob,
                const vector<vector<vector<double> > > &all_p) {
  // compute posterior probability at root node
  // (i.e. the init_state of all children)
  size_t lstate = all_paths[children[0]][site - 1].init_state;
  size_t rstate = all_paths[children[0]][site + 1].init_state;
  double p0 = root_trans_prob[lstate][0] * root_trans_prob[0][rstate];
  double p1 = root_trans_prob[lstate][1] * root_trans_prob[1][rstate];
  for (size_t idx = 0; idx < children.size(); ++idx) {
    p0 *= all_p[children[idx]][0][0];
    p1 *= all_p[children[idx]][0][1];
  }

  double root_p0 = p0 / (p0+p1);
  return root_p0;
}

void
downward_sampling(const vector<double> &triplet_rates,
                  const vector<size_t> &subtree_sizes,
                  const size_t site,
                  const vector<vector<Path> > &all_paths,
                  const vector<vector<double> > &root_trans_prob,
                  const vector<vector<vector<double> > > &all_p,
                  std::mt19937 &gen,
                  vector<Path> &new_path) {
  Path empty_path;
  new_path = vector<Path>(subtree_sizes.size(), empty_path);

  const size_t root_id = 0; //all_paths[0] is empty
  vector<size_t> children;
  get_children(root_id, subtree_sizes, children);

  // compute posterior probability at root node
  const double root_p0 = root_post_prob0(children, site, all_paths,
                                         root_trans_prob, all_p);

  // sample new root state
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  bool new_root_state = (unif(gen) > root_p0);
  for (size_t idx = 0; idx < children.size(); ++idx) {
    new_path[children[idx]].init_state = new_root_state;
    new_path[children[idx]].tot_time = all_paths[children[idx]][site - 1].tot_time;
  }

  // preorder traversal of the tree
  for (size_t node_id = 1; node_id < subtree_sizes.size(); ++node_id) {
    // collect relevant transition rates for each interval
    vector<vector<double> > interval_rates; //interval x transition type
    vector<double> interval_lengths;
    rates_on_branch(triplet_rates, all_paths[node_id][site-1],
                    all_paths[node_id][site+1],
                    interval_rates, interval_lengths);
    downward_sampling_branch(interval_rates, interval_lengths, subtree_sizes,
                             site, node_id, all_paths, all_p,
                             gen, new_path);
  }
}


////////////////////////////////////////////////////////////////////////////////
///////////////        Compute acceptance rate         /////////////////////////
////////////////////////////////////////////////////////////////////////////////

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


/*counterpart of downward_sampling_branch*/
void
proposal_prob_branch(const vector<vector<double> > &interval_rates,
                     const vector<double> &interval_lengths,
                     const vector<size_t> &subtree_sizes,
                     const size_t site,
                     const size_t node_id,
                     const vector<vector<Path> > &all_paths,
                     const vector<vector<vector<double> > > &all_p,
                     const Path &path,
                     double &prob) {

  const size_t n_intervals = interval_rates.size();
  const vector<double> jumps = path.jumps;

  size_t a = path.init_state; // state of one end
  double time_passed = 0;
  // end-conditioned sampling helpers
  vector<double> eigen_vals;
  vector<vector<double> > U;
  vector<vector<double> > Uinv;
  for (size_t m = 0; m < n_intervals - 1; ++m) {

    // compute conditional posterior probability
    vector<vector<double> > P; // transition matrix
    continuous_time_trans_prob_mat(interval_rates[m][0], interval_rates[m][1],
                                   interval_lengths[m], P);
    double p0 = (all_p[node_id][m+1][0] * P[a][0] /
                 all_p[node_id][m][a]);

    // check state at break point
    const double t = time_passed + interval_lengths[m];
    const size_t b = (size_t)(path.state_at_time(t)); // state of the other end
    prob *= (b == 0) ? p0: 1.0 - p0;

    // prepare helper values
    const CTMarkovModel ctmm(interval_rates[m]);

    // get jumps between time interval (time_passed, time_passed + interval_lengths[m])
    vector<double> subset_jumps;
    select_jumps(path.jumps, time_passed,
                 time_passed + interval_lengths[m], subset_jumps);

    // compute end_cond_sample prob
    prob *= end_cond_sample_prob(ctmm, a, b, interval_lengths[m], subset_jumps);

    // prepare for next interval
    time_passed += interval_lengths[m];
    a = b;
  }

  if (is_leaf(subtree_sizes[node_id])) {
    const size_t b = all_paths[node_id][site].end_state();
    const size_t m = n_intervals - 1;

    vector<double> subset_jumps;
    select_jumps(path.jumps, time_passed,
                 time_passed + interval_lengths[m], subset_jumps);

    // prepare helper values
    const CTMarkovModel ctmm(interval_rates[m]);

    // compute end_cond_sample prob
    prob *= end_cond_sample_prob(ctmm, a, b, interval_lengths[m], subset_jumps);
  }
  else {
    // last interval requires information from children
    vector<size_t> children;
    get_children(node_id, subtree_sizes, children);
    const size_t m = n_intervals - 1;

    // prepare helper values
    const CTMarkovModel ctmm(interval_rates[m]);

    // compute conditional posterior probability
    vector<vector<double> > P; // transition matrix
    ctmm.get_trans_prob_mat(interval_lengths[m], P);

    double p0 = 1.0;
    for (size_t idx = 0; idx < children.size(); ++idx)
      p0 *= all_p[children[idx]][0][0];
    p0 *= P[a][0]/all_p[node_id][m][a];

    size_t b = (size_t)(path.end_state());
    prob *= (b == 0) ? p0 : 1.0 - p0;

    // select jumps
    vector<double> subset_jumps;
    select_jumps(path.jumps, time_passed,
                 time_passed + interval_lengths[m], subset_jumps);

    // compute the end-conditioned sampling probability
    prob *= end_cond_sample_prob(ctmm, a, b, interval_lengths[m], subset_jumps);
  }
}


// compute proposal prob
void
proposal_prob(const vector<double> &triplet_rates,
              const vector<size_t> &subtree_sizes,
              const size_t site,
              const vector<vector<Path> > &all_paths,
              const vector<vector<double> > &root_trans_prob,
              const vector<vector<vector<double> > > &all_p,
              const vector<Path> &new_path,
              double &pro_old, double &pro_new) {

  // compute posterior probability at root node
  const size_t root_id = 0;
  vector<size_t> children;
  get_children(root_id, subtree_sizes, children);
  const double root_p0 = root_post_prob0(children, site, all_paths,
                                         root_trans_prob, all_p);

  pro_old = 1.0;
  pro_new = 1.0;

  pro_old *= all_paths[1][site].init_state ? 1.0 - root_p0 : root_p0;
  pro_new *= new_path[1].init_state ? 1.0 - root_p0 : root_p0;

  // preorder traversal of the tree
  for (size_t node_id = 1; node_id < subtree_sizes.size(); ++node_id) {
    // collect relevant transition rates for each interval
    vector<vector<double> > interval_rates;
    vector<double> interval_lengths;
    rates_on_branch(triplet_rates, all_paths[node_id][site - 1],
                    all_paths[node_id][site+1],
                    interval_rates, interval_lengths);
    proposal_prob_branch(interval_rates, interval_lengths,
                         subtree_sizes, site, node_id, all_paths,
                         all_p, new_path[node_id], pro_new);
    proposal_prob_branch(interval_rates, interval_lengths,
                         subtree_sizes, site, node_id, all_paths,
                         all_p, all_paths[node_id][site], pro_old);
  }

  //  cerr << pro_old << "\t" << pro_new << endl;

}


static double
log_lik_ratio(const vector<double> &rates,
              const PathContextStat &pcs_num,
              const PathContextStat &pcs_denom) {
  double result = 0.0;
  for (size_t i = 0; i < 8; ++i) {
    result += (pcs_num.jumps_in_context[i] -
               pcs_denom.jumps_in_context[i])*log(rates[i]) -
      (pcs_num.time_in_context[i] -
       pcs_denom.time_in_context[i]) * rates[i];
  }
  return result;
}





// compute acceptance rate
double
log_accept_rate(const EpiEvoModel &the_model, const TreeHelper &th,
                const size_t site,
                const vector<vector<Path> > &all_paths,
                const vector<vector<vector<double> > > &all_p,
                const vector<Path> &new_path) {
  double pro_old, pro_new;
  proposal_prob(the_model.triplet_rates, th.subtree_sizes,
                site, all_paths, the_model.init_T, all_p, new_path,
                pro_old, pro_new);

  double lr = log(pro_old) - log(pro_new);

  for (size_t i = 1; i < th.subtree_sizes.size(); ++i) {
    Path l = all_paths[i][site - 1];
    Path ll = all_paths[i][site - 2];
    Path r = all_paths[i][site + 1];
    Path rr = all_paths[i][site + 2];
    PathContextStat pcs_old(l, all_paths[i][site], r);
    PathContextStat pcs_new(l, new_path[i], r);
    PathContextStat pcs_old_l(ll, l, all_paths[i][site]);
    PathContextStat pcs_new_l(ll, l, new_path[i]);
    PathContextStat pcs_old_r(all_paths[i][site], r, rr);
    PathContextStat pcs_new_r(new_path[i], r, rr);

    lr += log_lik_ratio(the_model.triplet_rates, pcs_new, pcs_old) +
      log_lik_ratio(the_model.triplet_rates, pcs_new_l, pcs_old_l) +
      log_lik_ratio(the_model.triplet_rates, pcs_new_r, pcs_old_r);
  }

  return lr;
}


void
gibbs_site(const EpiEvoModel &the_model, const TreeHelper &th,
           const size_t site,
           vector<vector<Path> > &all_paths,
           std::mt19937 &gen,
           vector<Path> &new_path) {

  // collect relevant transition rates for each interval
  const size_t n_nodes = th.subtree_sizes.size();
  vector<vector<vector<double> > > all_interval_rates(n_nodes);
  vector<vector<double> > all_interval_lengths(n_nodes);
  for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
    rates_on_branch(the_model.triplet_rates,
                    all_paths[node_id][site - 1],
                    all_paths[node_id][site + 1],
                    all_interval_rates[node_id],
                    all_interval_lengths[node_id]);
  }

  // upward pruning and downward sampling
  vector<vector<vector<double> > > all_p;
  pruning(the_model.triplet_rates, th.subtree_sizes, site, all_paths,
          all_interval_rates, all_interval_lengths, all_p);

  downward_sampling(the_model.triplet_rates, th.subtree_sizes, site,
                    all_paths, the_model.init_T, all_p, gen, new_path);

  // acceptance rate
  const double log_acc_rate = log_accept_rate(the_model, th, site, all_paths,
                                              all_p, new_path);
  std::uniform_real_distribution<double> unif(0.0, 1.0);

  if (log_acc_rate >= 0 || unif(gen) < exp(log_acc_rate))
    for (size_t i = 1; i < th.subtree_sizes.size(); ++i)
      all_paths[i][site] = new_path[i];
}
