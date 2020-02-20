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

#include "ContinuousTimeMarkovModel.hpp"

#include <string>
#include <vector>
#include <ostream>
#include <sstream>
#include <cassert>
#include <cmath>

using std::string;
using std::vector;
using std::pair;
using std::endl;

static void
factor_rate_matrix(const double rate0, const double rate1,
                   vector<double> &eigen_vals,
                   two_by_two<double> &U,
                   two_by_two<double> &Uinv) {

  // Q = U*D*Uinv
  eigen_vals = vector<double>(2, 0.0);
  const double sum_rate = rate0 + rate1;
  eigen_vals[1] = -sum_rate;

  U(0, 0) = 1.0;
  U(0, 1) = rate0;
  U(1, 0) = 1.0;
  U(1, 1) = -rate1;

  Uinv(0, 0) = rate1 / sum_rate;
  Uinv(0, 1) = rate0 / sum_rate;
  Uinv(1, 0) = 1.0 / sum_rate;
  Uinv(1, 1) = -1.0 / sum_rate;
}


CTMarkovModel::CTMarkovModel(const vector<double> &rates) :
  rate0(rates.front()), rate1(rates.back()) {
  factor_rate_matrix(rate0, rate1, eigen_values, U, Uinv);
}

CTMarkovModel::CTMarkovModel(const pair<double, double> &rates) :
  rate0(rates.first), rate1(rates.second) {
  factor_rate_matrix(rate0, rate1, eigen_values, U, Uinv);
}

CTMarkovModel::CTMarkovModel(const double r0, const double r1) :
  rate0(r0), rate1(r1) {
  factor_rate_matrix(rate0, rate1, eigen_values, U, Uinv);
}


void
CTMarkovModel::get_trans_prob_mat(const double time_interval,
                                  two_by_two<double> &transition_matrix) const {

  assert(rate0 > 0 && rate1 > 0 && time_interval > 0);

  const double h = 1.0 / std::exp(time_interval * (rate0 + rate1));

  const double denominator = rate0 + rate1;
  transition_matrix(0, 0) = (rate0*h + rate1)/denominator;
  transition_matrix(0, 1) = 1.0 - transition_matrix(0, 0);

  transition_matrix(1, 1) = (rate0 + rate1*h)/denominator;
  transition_matrix(1, 0) = 1.0 - transition_matrix(1, 1);
}


string
CTMarkovModel::tostring() const {
  std::ostringstream oss;
  oss.precision(3);
  oss << "rates  =[" << rate0 << '\t' << rate1 << ']' << endl
      << "U      =[" << U(0, 0) << '\t' << U(0, 1) << ']' << endl
      << "        [" << U(1, 0) << '\t' << U(1, 1) << ']' << endl
      << "Uinv   =[" << Uinv(0, 0) << '\t' << Uinv(0, 1) << ']' << endl
      << "        [" << Uinv(1, 0) << '\t' << Uinv(1, 1) << ']' << endl
      << "lambda =[" << eigen_values[0] << '\t' << eigen_values[1] << ']';
  return oss.str();
}

std::ostream &
operator<<(std::ostream &os, const CTMarkovModel &ctmm) {
  return os << ctmm.tostring();
}


////////////////////////////////////////////////////////////////////////////////
///////////    SPECIALIZED CLASS "TwoStateCTMarkovModel" BELOW HERE    /////////
////////////////////////////////////////////////////////////////////////////////

double
TwoStateCTMarkovModel::get_trans_prob(const double time_interval,
                                      const size_t start_state,
                                      const size_t end_state) const {
  assert(rate0 > 0 && rate1 > 0 && time_interval > 0);
  const double h = exp(- time_interval * (rate0 + rate1));
  const double denominator = rate0 + rate1;
  const double prob =
    (start_state ? rate0 + rate1*h : rate0*h + rate1)/denominator;
  return (end_state == start_state ? prob : 1.0 - prob);
}

string
TwoStateCTMarkovModel::tostring() const {
  std::ostringstream oss;
  oss.precision(3);
  oss << "rates  =[" << rate0 << '\t' << rate1 << ']' << endl;
  return oss.str();
}

std::ostream &
operator<<(std::ostream &os, const TwoStateCTMarkovModel &tsctmm) {
  return os << tsctmm.tostring();
}


/* 2x2 rate matrix to transition prob. matrix */
void
continuous_time_trans_prob_mat(const double rate0, const double rate1,
                               const double time_interval,
                               two_by_two<double> &transition_matrix) {

  assert(rate0 > 0 && rate1 > 0 && time_interval > 0);

  const double h = 1.0/exp(time_interval*(rate0 + rate1));

  // ADS: this might be useful at some point here:
  // assert(h > 0.0);

  const double denominator = rate0 + rate1;
  transition_matrix(0, 0) = (rate0*h + rate1)/denominator;

  transition_matrix(0, 1) = 1.0 - transition_matrix(0, 0);

  transition_matrix(1, 1) = (rate0 + rate1*h)/denominator;
  transition_matrix(1, 0) = 1.0 - transition_matrix(1, 1);
}

////////////////////////////////////////////////////////////////////////////////
///////////////////       CTMM-RELATED STATISTICS     //////////////////////////
////////////////////////////////////////////////////////////////////////////////

void
expectation_J(const double r0, const double r1, const double T,
              two_by_two<double> &J0, two_by_two<double> &J1) {

  const double s = r0 + r1;
  const double p = r0 * r1;
  const double d = r1 - r0;
  const double e = exp(- s * T);

  const double C1 = d * (1 - e) / s;
  J0(0, 0) = p * ( T * (r1 - r0*e) - C1 ) / (s * (r1 + r0*e));
  J1(0, 0) = J0(0, 0);
  assert(J0(0, 0) >= 0);

  J0(1, 1) = p * ( T * (r0 - r1*e) + C1 ) / (s * (r0 + r1*e));
  J1(1, 1) = J0(1, 1);
  assert(J0(1, 1) >= 0);

  const double C2 = p * T * (1 + e) / (s * (1 - e));
  const double C3 = (r0 * r0 + r1 * r1) / (s * s);
  const double C4 = (2 * p) / (s * s);

  J0(0, 1) = C2 + C3;
  J1(0, 1) = C2 - C4;
  J0(1, 0) = J1(0, 1);
  J1(1, 0) = J0(0, 1);
  assert(J0(0, 1) >= 0);
  assert(J1(0, 1) >= 0);
}


void
expectation_D(const double r0, const double r1, const double T,
              two_by_two<double> &D0, two_by_two<double> &D1) {

  const double r00 = r0 * r0;
  const double r11 = r1 * r1;
  const double s = r0 + r1;
  const double p = r0 * r1;
  const double e = exp(- s * T);

  const double C1 = 2 * p * (1 - e) / s;
  D0(0, 0) = ((r11 + r00 * e) * T + C1) / (s * (r1 + r0 * e));
  D1(0, 0) = T - D0(0, 0);
  assert((D0(0, 0) >= 0) && (D0(0, 0) <= T));

  D1(1, 1) = ((r00 + r11 * e) * T + C1) / (s * (r0 + r1 * e));
  D0(1, 1) = T - D1(1, 1);
  assert((D1(1, 1) >= 0) && (D1(1, 1) <= T));

  const double C2 = (p - r00) * (1 - e) / s;
  D1(0, 1) = ((r00 - p * e) * T + C2) / (s * (r0 - r0 * e));
  D0(0, 1) = T - D1(0, 1);
  assert((D1(0, 1) >= 0) && (D1(0, 1) <= T));

  const double C3 = (p - r11) * (1 - e) / s;
  D0(1, 0) = ((r11 - p * e) * T + C3) / (s * (r1 - r1 * e));
  D1(1, 0) = T - D0(1, 0);
  assert((D0(1, 0) >= 0) && (D0(1, 0) <= T));
}
