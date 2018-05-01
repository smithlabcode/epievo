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

#include "EndCondSampler.hpp"
#include "StateSeq.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;

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

/* cdf function of end-conditioned time of first jump within the interval
   a: init state, b: end state, x: time of first jump
 */
double
cdf(const vector<double> &rates,
    const vector<double> &eigen_vals,
    const vector<vector<double> > &U,
    const vector<vector<double> > &Uinv,
    const vector<vector<double> > &PT,
    const double T, const size_t a, const size_t b, const double x) {
  const size_t a_bar = complement_state(a);
  const double a_rate = rates[a];
  double p = 0.0;
  for (size_t i = 0; i < 2; ++i) {
    const double scaler = U[a_bar][i] * Uinv[i][b];
    const double coeff = eigen_vals[i] + a_rate;
    const double expon = exp(T * eigen_vals[i])
                         - exp(T * eigen_vals[i] - x * coeff);
    p += scaler * (1.0 / coeff) * expon;
  }
  p *= a_rate/PT[a][b];
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
  const double w = line_search_cdf(rates, eigen_vals, U, Uinv, PT, T, a, b, rn);

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
      start_state = complement_state(start_state);
    }
    base_time += wait;
    tot -= wait;
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
