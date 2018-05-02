/* Copyright (C) 2018 University of Southern California
 *                    Jianghan Qu, Andrew D Smith and Xiaojing Ji
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
#include <algorithm>   // std::lower_bound,
#include <iostream>
#include <cmath>
#include <random>

#include "EndCondSampling.hpp"
#include "StateSeq.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;

static const double NUMERICAL_TOLERANCE = 1e-10;

// ADS: I think the MINWAIT value needs to be considered for
// elimination from the code.
static const double MINWAIT = 1e-8;

/* 2x2 rate matrix to transition prob. matrix */
void
trans_prob_mat(const double rate0, const double rate1,
               const double time_interval,
               vector<vector<double> > &transition_matrix) {

  assert(rate0 > 0 && rate1 > 0 && time_interval > 0);

  const double h = 1.0 / exp(time_interval * (rate0 + rate1));

  transition_matrix = vector<vector<double> >(2, vector<double>(2, 0.0));

  const double denominator = rate0 + rate1;
  transition_matrix[0][0] = (rate0*h + rate1)/denominator;
  transition_matrix[0][1] = 1.0 - transition_matrix[0][0];

  transition_matrix[1][1] = (rate0 + rate1*h)/denominator;
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
  const size_t a_bar = complement_state(a);
  for (size_t i = 0; i < 2; ++i)
    f += U[a_bar][i] * Uinv[i][b] * exp(T * eigen_vals[i]) * exp(-x * (eigen_vals[i] + rates[a]));

  f *= rates[a] / PT[a][b];
  return f;
}


/* This function is J_{aj} in Hobolth & Stone (2009), which appears on
   page 1209, is part of REMARK 4. In our case, with only two states,
   we don't need to handle the case of (lambda_j + Q_q = 0).
 */
static double
transformed_cdf_integrand_helper(const double total_time,  // T (HS2009)
                                 const double the_eig,     // lambda_j (HS2009)
                                 const double the_rate,    // Q_a (HS2009)
                                 const double proposed_time) { // t (HS2009)

  const double neg_eig_plus_rate = -1.0*(the_eig + the_rate);
  const double tot_time_mult_eig = total_time*the_eig;

  const double r =
    (exp(tot_time_mult_eig + proposed_time*neg_eig_plus_rate) -
     exp(tot_time_mult_eig))/neg_eig_plus_rate;

  assert(std::isfinite(r));

  return r;
}


/* This function computes the value in the summation of eqn (2.5) of
   Hobolth & Stone (2009). The coefficient of Q_{ai}/P_{ab}(T) is not
   included here, and the value "j" from that paper may only take
   values of 0 or 1 for our binary state space.
 */
double
transformed_cdf(const vector<double> &rates,
                const vector<double> &eigen_vals,
                const vector<vector<double> > &U,
                const vector<vector<double> > &Uinv,
                const double T, const size_t a, const size_t b, const double x) {

  const size_t a_bar = complement_state(a);
  const double a_rate = rates[a];

  const double integrand_0 =
    transformed_cdf_integrand_helper(T, eigen_vals[0], a_rate, x);
  const double intermediate_factor_0 = U[a_bar][0] * Uinv[0][b];

  const double integrand_1 =
    transformed_cdf_integrand_helper(T, eigen_vals[1], a_rate, x);
  const double intermediate_factor_1 = U[a_bar][1] * Uinv[1][b];

  return std::max(intermediate_factor_0*integrand_0 +
                  intermediate_factor_1*integrand_1,
                  std::numeric_limits<double>::min());
}


static double
cumulative_density(const vector<double> &rates,
                   const vector<double> &eigen_vals,
                   const vector<vector<double> > &U,
                   const vector<vector<double> > &Uinv,
                   const vector<vector<double> > &PT,
                   const double T, const size_t a, const size_t b,
                   const double x) {

  return transformed_cdf(rates, eigen_vals, U, Uinv, T, a, b, x)*
    rates[a]/PT[a][b];
}


static double
bisection_search_cumulative_density(const vector<double> &rates,
                                    const vector<double> &eigen_vals,
                                    const vector<vector<double> > &U,
                                    const vector<vector<double> > &Uinv,
                                    const vector<vector<double> > &PT,
                                    const double T,
                                    const size_t a, const size_t b,
                                    const double target) {
  double lo = 0.0;
  double lo_val = 0.0; // equals transformed_cdf for a value of 0.0

  double hi = T;
  double hi_val = transformed_cdf(rates, eigen_vals, U, Uinv, T, a, b, hi);

  // ADS: the target cdf value sampled from (0, x), where x is the
  // cumulative density for the -end-conditioned exponential is
  // transformed here to avoid the operations on mi_val below each
  // iteration, which also keeps values more centered.
  const double transformed_target = target/(rates[a]/PT[a][b]);

  assert(lo_val <= transformed_target && hi_val >= transformed_target);

  while (hi - lo > NUMERICAL_TOLERANCE) {
    const double mi = (lo + hi)/2.0;
    const double mi_val =
      transformed_cdf(rates, eigen_vals, U, Uinv, T, a, b, mi);
    if (mi_val >= transformed_target) {
      hi = mi;
      hi_val = mi_val;
    }
    else {
      lo = mi;
      lo_val = mi_val;
    }
  }
  return (lo + hi)/2.0;
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

  if (T < 2*std::numeric_limits<double>::min())
    return (a == b) ? T : T/2.0;

  vector<vector<double> > PT;  // PT = exp(QT)
  trans_prob_mat(rates[0], rates[1], T, PT);

  const double pr_no_jump = (a == b) ? exp(-rates[a]*T)/PT[a][a] : 0.0;

  std::uniform_real_distribution<double> unif(0.0, 1.0);
  if (a == b && unif(gen) < pr_no_jump)
    return T;

  // x~pdf <=> CDF(x)~Unif(0, 1)
  const double lowerbound = 0.0; // smallest possible cdf...
  // avoid calling cumulative_density for a value of 0.0

  const double upperbound =
    cumulative_density(rates, eigen_vals, U, Uinv, PT, T, a, b, T);

  std::uniform_real_distribution<double> unif_lu(lowerbound, upperbound);
  const double sampled_cdf_value = unif_lu(gen);

  const double w =
    bisection_search_cumulative_density(rates, eigen_vals, U, Uinv,
                                        PT, T, a, b, sampled_cdf_value);
  return w;
}


/* Endpoint-conditioned sampling of path witin time interval T */
void
end_cond_sample(const vector<double> rates,
                const vector<double> eigen_vals,
                const vector<vector<double> > U,
                const vector<vector<double> > Uinv,
                const size_t start_state, const size_t end_state, const double T,
                std::mt19937 &gen, vector<double> &jump_times) {

  jump_times.clear();

  size_t current_state = start_state;
  double consumed_time =
    end_cond_sample_first_jump(rates, eigen_vals, U, Uinv,
                               current_state, end_state, T, gen);

  // ADS: the use of NUMERICAL_TOLERANCE below should be checked. We
  // need to make sure that the sampling of a jump will return exactly
  // the total time interval when it should, and not some
  // approximation to it.
  while (T - consumed_time > NUMERICAL_TOLERANCE) {
    jump_times.push_back(consumed_time);
    current_state = complement_state(current_state);
    consumed_time +=
      end_cond_sample_first_jump(rates, eigen_vals, U, Uinv, current_state,
                                 end_state, T - consumed_time, gen);
  }
}


/* Endpoint-conditioned probability density*/
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

    a = complement_state(a);
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
