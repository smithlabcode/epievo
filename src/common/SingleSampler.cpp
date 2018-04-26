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

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;

static const double MINWAIT = 1e-8;

/* 2x2 rate matrix to transition prob. matrix */
static void
trans_prob_mat(const double rate0, const double rate1,
               const double interval,
               vector<vector<double> > &transition_matrix) {

  assert(rate0 > 0 && rate1 > 0 && interval > 0);
  double h = 1.0 / exp(interval * (rate0 + rate1));

  transition_matrix = vector<vector<double> >(2, vector<double>(2, 0.0));

  transition_matrix[0][0] = (rate0 * h + rate1) / (rate0 + rate1);
  transition_matrix[0][1] = 1.0 - transition_matrix[0][0];

  transition_matrix[1][1] = (rate0 + rate1 * h) / (rate0 + rate1);
  transition_matrix[1][0] = 1.0 - transition_matrix[1][1];
}

void
decompose(const vector<double> &rates, // rate0 and rate1
          vector<double> &eigen_vals,
          vector<vector<double> > &U,
          vector<vector<double> > &Uinv) {

  // Q = U*D*Uinv
  eigen_vals = vector<double>(2, 0.0);
  const double sum_rate = rates[0] + rates[1];
  eigen_vals[1] = -sum_rate;
  U = vector<vector<double> >(2, vector<double>(2, 0.0));
  U[0][0] = 1.0;
  U[0][1] = rates[0];
  U[1][0] = 1;
  U[1][1] = -rates[1];
  Uinv = vector<vector<double> >(2, vector<double>(2, 0.0));
  Uinv[0][0] = rates[1] / sum_rate;
  Uinv[0][1] = rates[0] / sum_rate;
  Uinv[1][0] = 1.0 / sum_rate;
  Uinv[1][1] = -1.0 / sum_rate;
}

/* pdf function of end-conditioned time of first jump within the interval*/
double
pdf(const vector<double> &rates,
    const vector<double> &eigen_vals,
    const vector<vector<double> > &U,
    const vector<vector<double> > &Uinv,
    const vector<vector<double> > &PT,
    const double T, const size_t a, const size_t b, const double x) {
  double f = 0.0;
  const size_t abar = 1 - a;
  for (size_t i = 0; i < 2; ++i)
    f += U[abar][i] * Uinv[i][b] * exp(T * eigen_vals[i]) * exp(-x * (eigen_vals[i] + rates[a]));

  f *= rates[a] / PT[a][b];
  return f;
}

double
cdf(const vector<double> &rates,
    const vector<double> &eigen_vals,
    const vector<vector<double> > &U,
    const vector<vector<double> > &Uinv,
    const vector<vector<double> > &PT,
    const double T, const size_t a, const size_t b, const double x) {
  double p = 0.0;
  const size_t abar = 1 - a;
  for (size_t i = 0; i < 2; ++i) {
    const double scaler = U[abar][i] * Uinv[i][b] * exp(T * eigen_vals[i]);
    const double coeff = -(eigen_vals[i] + rates[a]);
    p += scaler * (1.0 / coeff) * (exp(x * coeff) - 1.0);
  }
  p *= rates[a]/PT[a][b];
  return p;
}

double
line_search_cdf(const vector<double> &rates,
                const vector<double> &eigen_vals,
                const vector<vector<double> > &U,
                const vector<vector<double> > &Uinv,
                const vector<vector<double> > &PT,
                const double T, const size_t a, const size_t b,
                const double target) {
  double lo = MINWAIT; // PRECISION
  double hi = T-MINWAIT;
  double mi = 0.5 * (lo + hi);

  while (hi - lo > MINWAIT) {  // PRECISION
    double lo_val = cdf(rates, eigen_vals, U, Uinv, PT, T, a, b, lo);
    double mi_val = cdf(rates, eigen_vals, U, Uinv, PT, T, a, b, mi);
    double hi_val = cdf(rates, eigen_vals, U, Uinv, PT, T, a, b, hi);
    assert(lo_val <= target && hi_val >= target);
    if (mi_val > target) {
      hi = mi;
      mi = 0.5 * (lo + hi);
    } else if (mi_val < target) {
      lo = mi;
      mi = 0.5 * (lo + hi);
    } else {
      return mi;
    }
  }
  return mi;
}

/* Continuous time Markov chian with rate matrix Q.
   Return the first jump time within (0, T) or T if no jumps,
   given state at time 0 being a, and state at T being b.
*/
double
end_cond_sample_first_jump(const vector<double> rates,
                           const vector<double> eigen_vals,
                           const vector<vector<double> > U,
                           const vector<vector<double> > Uinv,
                           const size_t a, const size_t b,
                           const double T, std::mt19937 &gen) {

  if (a == b && T <= 2*MINWAIT) return T;
  if (a != b && T <= 2*MINWAIT) return T/2;

  vector<vector<double> > PT;  // PT = exp(QT)
  trans_prob_mat(rates[0], rates[1], T, PT);

  const double pr_no_jump  = (a == b) ? (exp(-rates[a] * T) / PT[a][a]) : 0.0;

  std::uniform_real_distribution<double> unif(0.0, 1.0);
  if (a == b && unif(gen) < pr_no_jump)
    return T;

  // x~pdf <=> CDF(x)~Unif(0, 1)

  const double upperbound = cdf(rates, eigen_vals, U, Uinv, PT, T, a, b, T-MINWAIT);
  const double lowerbound = cdf(rates, eigen_vals, U, Uinv, PT, T, a, b, MINWAIT);

  std::uniform_real_distribution<double> unif_lu(lowerbound, upperbound);
  const double rn = unif_lu(gen);

  // now do line search to find x s.t. CDF(x) = rn
  double w = line_search_cdf(rates, eigen_vals, U, Uinv, PT, T, a, b, rn);
  return w;
}


/* Endpoint-conditioned sampling of path witin time interval T*/
void
end_cond_sample(const vector<double> rates,
                const vector<double> eigen_vals,
                const vector<vector<double> > U,
                const vector<vector<double> > Uinv,
                const size_t a, const size_t b, const double T,
                std::mt19937 &gen, vector<double> &jump_times) {

  jump_times.clear();
  double tot = T;
  size_t start_state = a;
  double base_time = 0;
  while (tot > 0) {
    double wait = end_cond_sample_first_jump(rates, eigen_vals, U, Uinv,
                                             start_state, b, tot, gen);

    if (wait < tot) {
      jump_times.push_back(base_time + wait);
      start_state = 1 - start_state;
    }
    base_time += wait;
    tot -= wait;
  }
}



////////////////////////////////////////////////////////////////////////////////
///////////////////           UPWARD PRUNING           /////////////////////////
////////////////////////////////////////////////////////////////////////////////

// rates are two muation rates for pattern x0y, and x1y
void
upward(const vector<double> &rates,
       const double T,
       const vector<double> &q, // auxiliary values
       vector<double> &p) {

  assert(T > 0);

  p = vector<double>(2, 0.0);
  vector<vector<double> > P; // transition matrix
  trans_prob_mat(rates[0], rates[1], T, P);
  for (size_t j = 0; j < 2; ++j) {
    for (size_t k = 0; k < 2; ++k) {
      p[j] += P[j][k]*q[k];
    }
  }
}


// collect rates and interval lengths
void
rates_on_branch(const vector<double> &triplet_rates,
                const Path &l, const Path &r,
                vector<vector<double> > &interval_rates,
                vector<double> &interval_lengths) {
  Environment env(l, r);
  const size_t n_intervals = env.left.size();
  interval_rates = vector<vector<double> > (n_intervals, vector<double>(2, 0.0));
  interval_lengths = vector<double>(n_intervals, 0.0);

  for (size_t i = 0; i < n_intervals; ++i) {
    const size_t pattern0 = 4 * (size_t)(env.left[i]) + (size_t)(env.right[i]);
    const size_t pattern1 = pattern0 + 2;
    interval_rates[i][0] = triplet_rates[pattern0];
    interval_rates[i][1] = triplet_rates[pattern1];
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
               const vector<vector<vector<double> > > &all_interval_rates,
               const vector<vector<double> > & all_interval_lengths,
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
          for (size_t k = 0; k < 2; ++k) {
            q = vector<double>(2, 1.0);
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

  size_t par_state = new_path[node_id].init_state;
  double time_passed = 0;
  // end-conditioned sampling helpers
  vector<double> eigen_vals;
  vector<vector<double> > U;
  vector<vector<double> > Uinv;
  for (size_t m = 0; m < n_intervals - 1; ++m) {

    // compute conditional posterior probability
    vector<vector<double> > P; // transition prob matrix
    trans_prob_mat(interval_rates[m][0], interval_rates[m][1],
                   interval_lengths[m], P);
    double p0 = (all_p[node_id][m+1][0] * P[par_state][0] /
                 all_p[node_id][m][par_state]);

    // generate random state at break point
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    bool new_state = (unif(gen) > p0);

    // prepare helper values
    decompose(interval_rates[m], eigen_vals, U, Uinv);

    // generate path
    vector<double> jump_times;
    end_cond_sample(interval_rates[m], eigen_vals, U, Uinv,
                    par_state, (size_t)(new_state),
                    interval_lengths[m], gen, jump_times);

    // append jump_times to new_path
    for (size_t i = 0; i < jump_times.size(); ++i) {
      new_path[node_id].jumps.push_back(time_passed + jump_times[i]) ;
    }

    // prepare for next interval
    time_passed += interval_lengths[m];
    par_state = new_state;
  }

  if (is_leaf(subtree_sizes[node_id])) {
    const size_t leaf_state = all_paths[node_id][site].end_state();
    // generate path
    vector<double> jump_times;
    const size_t m = n_intervals - 1;
    // prepare helper values
    decompose(interval_rates[m], eigen_vals, U, Uinv);
    end_cond_sample(interval_rates[m], eigen_vals, U, Uinv, par_state, leaf_state,
                    interval_lengths[m], gen, jump_times);
    // append jump_times to new_path
    for (size_t i = 0; i < jump_times.size(); ++i) {
      new_path[node_id].jumps.push_back(time_passed + jump_times[i]) ;
    }

    assert(new_path[node_id].end_state() == leaf_state);

  } else {
    // last interval requires information from children
    vector<size_t> children;
    get_children(node_id, subtree_sizes, children);
    const size_t m = n_intervals - 1;
    // compute conditional posterior probability
    vector<vector<double> > P; // transition matrix
    trans_prob_mat(interval_rates[m][0], interval_rates[m][1],
                   interval_lengths[m], P);
    double p0 = 1.0;
    for (size_t idx = 0; idx < children.size(); ++idx)
      p0 *= all_p[children[idx]][0][0];
    p0 *= P[par_state][0] / all_p[node_id][m][par_state];

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
    decompose(interval_rates[m], eigen_vals, U, Uinv);
    // generate path
    vector<double> jump_times;
    end_cond_sample(interval_rates[m], eigen_vals, U, Uinv,
                    par_state, (size_t)(new_state),
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

double
end_cond_sample_prob(const vector<double> rates,
                     const vector<double> eigen_vals,
                     const vector<vector<double> > U,
                     const vector<vector<double> > Uinv,
                     size_t a, const size_t b, double T,
                     vector<double> jump_times) {
  double p = 1.0;

  // jump_times are between 0 and T
  assert(jump_times.size() == 0 || jump_times.back() < T);
  assert(jump_times.size() > 0 || a == b);

  while (jump_times.size()) {
     vector<vector<double> > PT;  // PT = exp(QT)
     trans_prob_mat(rates[0], rates[1], T, PT);
     p *= pdf(rates, eigen_vals, U, Uinv, PT, T, a, b, jump_times[0]);

    a = 1 - a;
    const double w = jump_times[0];
    T = T - w;
    jump_times.erase(jump_times.begin());
    for (size_t i = 0; i < jump_times.size(); ++i)
      jump_times[i] -= w;
  }

  assert(a == b);
  vector<vector<double> > PT;  // PT = exp(QT)
  trans_prob_mat(rates[0], rates[1], T, PT);
  const double pr_no_jump = exp(-rates[a] * T) / PT[a][a];
  return p * pr_no_jump;
}


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
    trans_prob_mat(interval_rates[m][0], interval_rates[m][1],
                   interval_lengths[m], P);
    double p0 = (all_p[node_id][m+1][0] * P[a][0] /
                 all_p[node_id][m][a]);

    // check state at break point
    const double t = time_passed + interval_lengths[m];
    const size_t b = (size_t)(path.state_at_time(t)); // state of the other end
    prob *= (b == 0) ? p0: 1.0 - p0;

    // prepare helper values
    decompose(interval_rates[m], eigen_vals, U, Uinv);

    // get jumps between time interval (time_passed, time_passed + interval_lengths[m])
    vector<double> subset_jumps;
    select_jumps(path.jumps, time_passed,
                 time_passed + interval_lengths[m], subset_jumps);

    // compute end_cond_sample prob
    prob *= end_cond_sample_prob(interval_rates[m], eigen_vals, U, Uinv,
                                 a, b, interval_lengths[m], subset_jumps);

    // prepare for next interval
    time_passed += interval_lengths[m];
    a = b;
  }

  if (is_leaf(subtree_sizes[node_id])) {
    const size_t b = all_paths[node_id][site].end_state();
    const size_t m = n_intervals - 1;
    // prepare helper values
    decompose(interval_rates[m], eigen_vals, U, Uinv);

    vector<double> subset_jumps;
    select_jumps(path.jumps, time_passed,
                 time_passed + interval_lengths[m], subset_jumps);
    // compute end_cond_sample prob
    prob *= end_cond_sample_prob(interval_rates[m], eigen_vals, U, Uinv,
                                 a, b, interval_lengths[m], subset_jumps);
  } else {
    // last interval requires information from children
    vector<size_t> children;
    get_children(node_id, subtree_sizes, children);
    const size_t m = n_intervals - 1;

    // compute conditional posterior probability
    vector<vector<double> > P; // transition matrix
    trans_prob_mat(interval_rates[m][0], interval_rates[m][1],
                   interval_lengths[m], P);
    double p0 = 1.0;
    for (size_t idx = 0; idx < children.size(); ++idx)
      p0 *= all_p[children[idx]][0][0];
    p0 *= P[a][0]/all_p[node_id][m][a];

    size_t b = (size_t)(path.end_state());
    prob *= (b == 0) ? p0 : 1.0 - p0;

    // prepare helper values
    decompose(interval_rates[m], eigen_vals, U, Uinv);

    // select jumps
    vector<double> subset_jumps;
    select_jumps(path.jumps, time_passed,
                 time_passed + interval_lengths[m], subset_jumps);

    // compute end_cond_sample prob
    prob *= end_cond_sample_prob(interval_rates[m], eigen_vals, U, Uinv,
                                 a, b, interval_lengths[m], subset_jumps);
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
log_accept_rate(const EpiEvoModel &the_model,
                const size_t site,
                const vector<vector<Path> > &all_paths,
                const vector<vector<vector<double> > > &all_p,
                const vector<Path> &new_path) {
  double pro_old, pro_new;
  proposal_prob(the_model.triplet_rates, the_model.subtree_sizes,
                site, all_paths, the_model.init_T, all_p, new_path,
                pro_old, pro_new);

  double lr = log(pro_old) - log(pro_new);

  for (size_t i = 1; i < the_model.subtree_sizes.size(); ++i) {
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
gibbs_site(const EpiEvoModel &the_model,
           const size_t site,
           vector<vector<Path> > &all_paths,
           std::mt19937 &gen,
           vector<Path> &new_path) {

  // collect relevant transition rates for each interval
  const size_t n_nodes = the_model.subtree_sizes.size();
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
  pruning(the_model.triplet_rates, the_model.subtree_sizes, site, all_paths,
          all_interval_rates, all_interval_lengths, all_p);

  downward_sampling(the_model.triplet_rates, the_model.subtree_sizes, site,
                    all_paths, the_model.init_T, all_p, gen, new_path);

  // acceptance rate
  double log_acc_rate = log_accept_rate(the_model, site, all_paths,
                                        all_p, new_path);
  std::uniform_real_distribution<double> unif(0.0, 1.0);

  if (log_acc_rate >= 0 || unif(gen) < exp(log_acc_rate))
    for (size_t i = 1; i < the_model.subtree_sizes.size(); ++i)
      all_paths[i][site] = new_path[i];


}
