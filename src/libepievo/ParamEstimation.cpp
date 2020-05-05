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
#include <array>
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <functional>

#include "PhyloTreePreorder.hpp"
#include "Path.hpp"
#include "EpiEvoModel.hpp"
#include "TreeHelper.hpp"
#include "TripletSampler.hpp" /* forward simulation */

#include "epievo_utils.hpp"

using std::vector;
using std::array;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::transform;
using std::begin;
using std::end;

template <class T, class U> static void
scale_mult(const T scale_factor, U &to_scale) {
  transform(begin(to_scale), end(to_scale), begin(to_scale),
            [&](const double x) {return x*scale_factor;});
}


template <class T, class U> static void
scale_div(const T scale_factor, U &to_scale) {
  transform(begin(to_scale), end(to_scale), begin(to_scale),
            [&](const double x) {return x/scale_factor;});
}


static void
add_sufficient_statistics(const vector<Path> &paths,
                          vector<double> &J, vector<double> &D) {
  // iterate over sites with valid triples (i.e. not the first and last)
  const size_t n_sites = paths.size();
  assert(n_sites > 0); // make sure the root isn't processed
  for (size_t i = 1; i < n_sites - 1; ++i)
    add_sufficient_statistics(paths[i-1], paths[i], paths[i+1], J, D);
}


void
get_sufficient_statistics(const vector<vector<Path> > &all_paths,
                          vector<double> &J, vector<double> &D) {
  static const size_t n_triplets = 8;
  // below result but faster than J = vector<double>(n_triplets, 0.0)
  J.clear();
  J.resize(n_triplets, 0.0);
  D.clear();
  D.resize(n_triplets, 0.0);

  // iterate over nodes, starting at 1 to avoid root
  for (size_t i = 1; i < all_paths.size(); ++i)
    add_sufficient_statistics(all_paths[i], J, D);
}

/* branches in likelihood and gradient ascent steps are scalers that
   should be multiplied to input branch lengths after the optimization
   is over */
void
get_sufficient_statistics(const vector<vector<Path> > &all_paths,
                          vector<vector<double> > &J,
                          vector<vector<double> > &D) {

  static const size_t n_triplets = 8;
  const size_t n_branches = all_paths[0].size();
  J.resize(n_branches);
  D.resize(n_branches);
  for (size_t i = 0; i < n_branches; ++i) {
    J[i].clear();
    J[i].resize(n_triplets, 0.0);
    D[i].clear();
    D[i].resize(n_triplets, 0.0);
  }

  // iterate over nodes, starting at 1 to avoid root
  for (size_t i = 1; i < all_paths.size() - 1; ++i)
    for (size_t b = 1; b < n_branches; ++b)
      add_sufficient_statistics(all_paths[i-1][b],
                                all_paths[i][b],
                                all_paths[i+1][b], J[b], D[b]);
}

void
get_root_frequencies(const vector<vector<Path>> &all_paths,
                     two_by_two &counts) {
  assert(!all_paths.empty() && all_paths[0].size() > 1);

  counts.reset();
  size_t prev = all_paths[0][1].init_state;
  for (size_t i = 1; i < all_paths.size(); ++i) {
    const size_t curr = all_paths[i][1].init_state;
    counts(prev, curr)++;
    prev = curr;
  }
}


static double
log_likelihood(const vector<double> &J, const vector<double> &D,
               const array<double, 8> &rates) {

  static const size_t n_triplets = 8;
  assert(J.size() == n_triplets && D.size() == n_triplets);

  double ll = 0;
  for (size_t i = 0; i < n_triplets; ++i)
    ll += J[i]*log(rates[i]) - D[i]*rates[i];

  return ll;
}


/* compute gradients wrt log(rate[i]) */
static void
add_to_gradient(const vector<double> &J, const vector<double> &D,
                const array<double, 8> &rates, vector<double> &gradient) {
  static const size_t n_params = 8;
  assert(gradient.size() == n_params);

  /* ordering of parameters:
     000 -> 010: 0 birth
     010 -> 000: 2 death
     001 -> 011: 1 expansion
     011 -> 001: 3 contraction
     100 -> 110: 4 expansion
     110 -> 100: 6 contraction
     101 -> 111: 5 merging
     111 -> 101: 7 splitting (not a free parameter)
  */
  const double factor_111 = (J[7] - D[7]*rates[7]);

  // [0] BIRTH: parameter 0 corresponds to 000->010
  gradient[0] += J[0] - D[0]*rates[0] + factor_111;

  // [2] DEATH: parameter 2 corresponds to 010->000
  gradient[2] += J[2] - D[2]*rates[2] - factor_111;

  // [1, 4] EXPANSION: parameter 1 corresponds to 001->011
  gradient[1] += J[1] + J[4] - (D[1] + D[4])*rates[1] - 2*factor_111;
  gradient[4] = gradient[1]; // expansion in other direction: 100->110

  // [3, 6] CONTRACTION: parameter 3 corresponds to 011->001
  gradient[3] += J[3] + J[6] - (D[3] + D[6])*rates[3] + 2*factor_111;
  gradient[6] = gradient[3]; // contraction in other direction: 110->100

  // [5] MERGING: parameter 2 corresponds to 101->111
  gradient[5] += J[5] - D[5]*rates[5] + factor_111;

  // [7] SPLITTING: parameter 7 corresponds to 111->101
  /* gradient[7] remains 0, since rates[7] is determined by other rates */
}

/* compute gradients wrt log(rate[i]) */
static void
get_gradient(const vector<double> &J, const vector<double> &D,
             const array<double, 8> &rates, vector<double> &gradient) {
  static const size_t n_params = 8;
  gradient.clear();
  gradient.resize(n_params, 0.0);
  add_to_gradient(J, D, rates, gradient);
}


/* Propose a set of new candidate rate by using the given gradient
 * vector and step size.
 */
static void
candidate_rates(const double step_size,
                const vector<double> &gradient,
                const array<double, 8> &rates,
                array<double, 8> &updated_rates) {

  static const size_t n_rates = 8;

  for (size_t i = 0; i < n_rates - 1; ++i)
    updated_rates[i] = exp(log(rates[i]) + gradient[i]*step_size);

  // final rate is in terms of other rates
  updated_rates[n_rates - 1] =
    exp(/**/log(updated_rates[0])   // 000 -> 010 (once, numerator)
        +   log(updated_rates[5])   // 101 -> 111 (once, numerator)
        + 2*log(updated_rates[3])   // 011 -> 001 (twice,numerator)
        -   log(updated_rates[2])   // 010 -> 000 (once, denominator)
        - 2*log(updated_rates[1])); // 001 -> 011 (twice,denominator)
}


/* This function calls the candidate_rates function to obtain new
   rates, and then updates the branch lengths correspondingly.
*/
static void
candidate_branches(const vector<vector<double> > &J,
                   const vector<vector<double> > &D,
                   const array<double, 8> &rates,
                   vector<double> &updated_branches) {

  assert(J.size() == D.size());
  const size_t n_branches = D.size();

  updated_branches.resize(n_branches);
  for (size_t b = 1; b < n_branches; ++b) {
    const double num = accumulate(begin(J[b]), end(J[b]), 0.0);
    const double denom = inner_product(begin(D[b]), end(D[b]),
                                       begin(rates), 0.0);
    updated_branches[b] = num/denom;
  }
}


static double
get_starting_step_size(const vector<double> &gradient) {
  double l1_norm = 0.0;
  for (auto &&i : gradient)
    l1_norm += fabs(i);
  return 1.0/l1_norm; // ADS: is this magic?
}


/* makes a single step of gradient ascent; identifies a new set of
 * rates using the summary statistics in J and D. The return value
 * indicates if the ascent was able to find a better likelihood.
 */
static bool
gradient_ascent(const double param_tol,
                const vector<double> &J, const vector<double> &D,
                const double llh, const array<double, 8> &rates,
                double &updated_llh, array<double, 8> &updated_rates) {

  assert(llh == log_likelihood(J, D, rates));

  /* compute llh and gradient */
  vector<double> gradient;
  get_gradient(J, D, rates, gradient);

  double step_size = get_starting_step_size(gradient);
  updated_llh = std::numeric_limits<double>::lowest();
  while (updated_llh < llh && step_size > param_tol) {
    candidate_rates(step_size, gradient, rates, updated_rates);
    updated_llh = log_likelihood(J, D, updated_rates);
    step_size *= 0.5;
  }
  return (updated_llh > llh);
}


static double
estimate_rates(const double param_tol,
               const vector<double> &J, const vector<double> &D,
               const array<double, 8> &input_rates,
               array<double, 8> &rates) {

  double llh = log_likelihood(J, D, input_rates);
  std::copy(begin(input_rates), end(input_rates), begin(rates));

  array<double, 8> tmp_rates;
  std::copy(begin(rates), end(rates), begin(tmp_rates));
  double tmp_llh = llh;
  while (gradient_ascent(param_tol, J, D, llh, rates, tmp_llh, tmp_rates)) {
    llh = tmp_llh;
    rates.swap(tmp_rates);
  }
  return llh;
}


static double
estimate_rates(const double param_tol,
               const vector<vector<double> > &J,
               const vector<vector<double> > &D,
               const array<double, 8> &input_rates,
               array<double, 8> &rates) {

  const size_t n_rates = J.back().size();
  vector<double> J_clps(n_rates, 0.0), D_clps(n_rates, 0.0);
  for (size_t b = 1; b < J.size(); ++b) {
    for (size_t i = 0; i < n_rates; ++i) {
      J_clps[i] += J[b][i];
      D_clps[i] += D[b][i];
    }
  }
  return estimate_rates(param_tol, J_clps, D_clps, input_rates, rates);
}


void
set_one_change_per_site_per_unit_time(array<double, 8> &rates,
                                      vector<double> &branches) {

  /* scale rates and branches to have unit branch length corresponding
     to one expected transition per site */
  const double the_rate_factor = rate_scaling_factor(rates);

  // scale up the branch lengths
  scale_mult(the_rate_factor, branches);

  // scale up the D values (total time spent in each state)
  // scale_mult(the_rate_factor, D);

  // scale down the rates
  scale_div(the_rate_factor, rates);
}


double
estimate_rates(const bool VERBOSE, const double param_tol,
               const vector<vector<double> > &J,
               const vector<vector<double> > &D,
               EpiEvoModel &the_model) {

  if (VERBOSE)
    cerr << "[ESTIMATING PARAMETERS]" << endl;

  array<double, 8> updated_rates;
  const double llh =
    estimate_rates(param_tol, J, D, the_model.triplet_rates, updated_rates);

  the_model.rebuild_from_triplet_rates(updated_rates);

  return llh;
}


double
estimate_rates(const bool VERBOSE, const double param_tol,
               const vector<vector<Path> > &all_paths,
               EpiEvoModel &the_model) {
  if (VERBOSE)
    cerr << "[COMPUTING INITIAL SUFFICIENT STATISTICS]" << endl;
  vector<vector<double>> J, D;
  get_sufficient_statistics(all_paths, J, D);

  return estimate_rates(VERBOSE, param_tol, J, D, the_model);
}


void
scale_jump_times(vector<vector<Path> > &all_paths, const TreeHelper &th) {
  for (size_t i = 0; i < all_paths.size(); ++i) {
    for (size_t b = 1; b < all_paths[i].size(); ++b) {
      const double scale = th.branches[b] / all_paths[i][b].tot_time;
      for (size_t j = 0; j < all_paths[i][b].jumps.size(); ++j) {
        all_paths[i][b].jumps[j] *= scale;
      }
      all_paths[i][b].tot_time = th.branches[b];
    }
  }
}


double
estimate_rates_and_branches(const bool VERBOSE, const double param_tol,
                            const vector<vector<double> > &J,
                            const vector<vector<double> > &D,
                            TreeHelper &th, EpiEvoModel &the_model) {

  if (VERBOSE)
    cerr << "[ESTIMATING PARAMETERS AND BRANCHES]" << endl;

  // stage 1: update rates
  array<double, 8> updated_rates;
  estimate_rates(param_tol, J, D, the_model.triplet_rates, updated_rates);

  // stage 2: update branches
  vector<double> updated_branches(th.branches);
  vector<double> branch_scale(th.branches.size(), 1.0);
  candidate_branches(J, D, updated_rates, branch_scale);
  transform(begin(branch_scale), end(branch_scale),
            begin(th.branches), begin(updated_branches),
            std::multiplies<double>());

  set_one_change_per_site_per_unit_time(updated_rates, updated_branches);
  the_model.rebuild_from_triplet_rates(updated_rates);

  // update branch lenths
  th.branches = updated_branches;

  // calculate new llh
  static const size_t n_triplets = 8;
  vector<double> J_collapsed(n_triplets, 0.0);
  vector<double> D_collapsed(n_triplets, 0.0);

  for (size_t b = 1; b < th.n_nodes; b++)
    for (size_t i = 0; i < D[b].size(); i++) {
      J_collapsed[i] += J[b][i];
      D_collapsed[i] += branch_scale[b]*D[b][i];
    }

  return log_likelihood(J_collapsed, D_collapsed, updated_rates);
}


double
estimate_rates_and_branches(const bool VERBOSE, const double param_tol,
                            vector<vector<Path>> &all_paths,
                            TreeHelper &th,
                            EpiEvoModel &the_model) {

  if (VERBOSE)
    cerr << "[NORMALIZING ALL PATH LENGTHS]" << endl;
  for (size_t b = 1; b < all_paths.size(); ++b)
    for (size_t i = 0; i < all_paths[b].size(); ++i) {
      all_paths[b][i].scale_to_unit_length(all_paths[b][i].tot_time);
      all_paths[b][i].tot_time = 1.0;
    }

  // get initial values of sufficient statistics
  if (VERBOSE)
    cerr << "[COMPUTING INITIAL SUFFICIENT STATISTICS]" << endl;
  vector<vector<double>> J, D;
  get_sufficient_statistics(all_paths, J, D);

  return estimate_rates_and_branches(VERBOSE, param_tol, J, D, th, the_model);
}
