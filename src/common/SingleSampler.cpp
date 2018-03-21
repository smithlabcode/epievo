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
#include <algorithm>
#include <iostream>
#include <cmath>
#include <random>

#include "SingleSampler.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;


/* 2x2 rate matrix to transition prob. matrix */
static void
trans_prob_mat(double rate0, double rate1,
               double interval,
               vector<vector<double> > &transition_matrix) {
  assert(rate0 > 0 && rate1 > 0 && interval > 0);
  double h = 1.0/exp(interval*(rate0 + rate1));

  transition_matrix = vector<vector<double> >(2, vector<double>(2,0.0));

  transition_matrix[0][0] = (rate0*h + rate1)/(rate0 + rate1);
  transition_matrix[0][1] = 1.0 - transition_matrix[0][0];

  transition_matrix[1][1] = (rate0 + rate1*h)/(rate0 + rate1);
  transition_matrix[1][0] = 1.0 - transition_matrix[1][1];
}



inline double
pdf(const double coeff0, const double coeff1,
    const double exp_coeff0, const double exp_coeff1,
    const double scaler, const double x) {

  return (coeff0 / exp_coeff0 * (exp(exp_coeff0 * x) - 1) +
          coeff1 / exp_coeff1 * (exp(exp_coeff1 * x) - 1))/scaler;
};




double
sample_from_jump_time_pdf(const vector<vector<double> > &Q,
                          const size_t a, const size_t b,
                          const double T, std::mt19937 &gen) {
  double coeff0 = 0;
  double coeff1 = 0;
  double exp_coeff0 = 0;
  double exp_coeff1 = 0;

  /* the pdf is: 1/scaler*[coeff0 * exp(exp_coeff0 * t) +  coeff1 * exp(exp_coeff1 * t)] */
  const double q0 = Q[0][1];
  const double q1 = Q[1][0];
  if (a == 0 && b == 0) {
    double denom = q1 + q0 * exp( - (q0 + q1) * T);
    coeff0 = pow(q0, 2) / denom;
    coeff1 = (-pow(q0, 2)) * exp(- (q0 + q1) * T) / denom;
    exp_coeff0 = -q0;
    exp_coeff1 = q1;
  }

  if (a == 0 && b == 1) {
    double denom = 1.0 - exp(- (q0 + q1) * T);
    coeff0 = q0/denom;
    coeff1 = q1 * exp(- (q0 + q1) * T) / denom;
    exp_coeff0 = -q0;
    exp_coeff1 = q1;
  }

  if (a == 1 && b == 0) {
    double denom = 1.0 - exp(- (q0 + q1) * T);
    coeff0 = q1 / denom;
    coeff1 = q0 * exp(- (q0 + q1) * T) / denom;
    exp_coeff0 = -q1;
    exp_coeff1 = q0;
  }

  if (a == 1 && b == 1) {
    double denom = q0 + q1 * exp( - (q0 + q1) * T);
    coeff0 = pow(q1, 2) / denom;
    coeff1 = (-pow(q1, 2)) * exp(- (q0 + q1) * T) / denom;
    exp_coeff0 = -q1;
    exp_coeff1 = q0;
  }

  const double scaler = coeff0 / exp_coeff0 * (exp(exp_coeff0 * T) - 1) +
    coeff1 / exp_coeff1 * (exp(exp_coeff1 * T) - 1);

  /* x~pdf <=> CDF(x)~Unif(0,1)*/
  std::uniform_real_distribution<double> unif(0.0,1.0);
  double rn = unif(gen);
  // now do line search to find x s.t. CDF(x) = rn
  double lo = 0;
  double hi = T;
  double mi = 0.5 * (lo + hi);
  while (hi-lo > 1e-5) {
    double lo_val = pdf(coeff0, coeff1, exp_coeff0, exp_coeff1, scaler, lo);
    double mi_val = pdf(coeff0, coeff1, exp_coeff0, exp_coeff1, scaler, mi);
    double hi_val = pdf(coeff0, coeff1, exp_coeff0, exp_coeff1, scaler, hi);
    assert (lo_val < rn && hi_val > rn);
    if (mi_val > rn) {
      hi = mi;
      mi = 0.5 * (lo + hi);
    } else if (mi_val < rn) {
      lo = mi;
      mi = 0.5 * (lo + hi);
    } else {
      return mi;
    }
  }
  return mi;
}



/* Continuous time Markov chian with rate matrix Q.
   Return the first jump time within (0,T) or T if no jumps,
   given state at time 0 being a, and state at T being b.
*/
double
end_cond_sample_first_jump(const vector<vector<double> > &Q,
                           const size_t a, const size_t b,
                           const double T, std::mt19937 &gen) {
  assert (a <= 1 && b <= 1);

  vector<vector<double> > PT;
  const double rate0 = Q[0][1];
  const double rate1 = Q[1][0];
  trans_prob_mat(rate0, rate1, T, PT);

  std::uniform_real_distribution<double> unif(0.0,1.0);
  double w = 0;
  if (a == b) {
    double pr_no_jump = exp(Q[a][a]*T - log(PT[a][a]));
    if (unif(gen) < pr_no_jump)
      return T;
  }

  w = sample_from_jump_time_pdf(Q, a, b, T, gen);
  cerr << "sampling w=" << w << endl;
  return w;
}

/* Endpoint-conditioned sampling of path witin time interval T*/
void
end_cond_sample(const std::vector<std::vector<double> > &Q,
                const size_t a, const size_t b, const double T,
                std::mt19937 &gen, vector<double> &jump_times) {
  jump_times = vector<double>(1, 0.0);
  double tot = T;
  size_t start_state = a;
  while (tot > 0) {
    double wait = end_cond_sample_first_jump(Q, start_state, b, tot, gen);
    assert(wait > 0);
    jump_times.push_back(jump_times.back() + wait);
    start_state = 1 - start_state;
    tot -= wait;
  }
}
