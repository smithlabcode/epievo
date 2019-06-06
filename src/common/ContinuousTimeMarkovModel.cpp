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
                   vector<vector<double> > &U,
                   vector<vector<double> > &Uinv) {

  // Q = U*D*Uinv
  eigen_vals = vector<double>(2, 0.0);
  const double sum_rate = rate0 + rate1;
  eigen_vals[1] = -sum_rate;

  U = two_by_two(2, vector<double>(2, 0.0));
  U[0][0] = 1.0;
  U[0][1] = rate0;
  U[1][0] = 1;
  U[1][1] = -rate1;

  Uinv = two_by_two(2, vector<double>(2, 0.0));
  Uinv[0][0] = rate1 / sum_rate;
  Uinv[0][1] = rate0 / sum_rate;
  Uinv[1][0] = 1.0 / sum_rate;
  Uinv[1][1] = -1.0 / sum_rate;
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
                                  two_by_two &transition_matrix) const {

  assert(rate0 > 0 && rate1 > 0 && time_interval > 0);

  const double h = 1.0 / std::exp(time_interval * (rate0 + rate1));

  transition_matrix = vector<vector<double> >(2, vector<double>(2, 0.0));

  const double denominator = rate0 + rate1;
  transition_matrix[0][0] = (rate0*h + rate1)/denominator;
  transition_matrix[0][1] = 1.0 - transition_matrix[0][0];

  transition_matrix[1][1] = (rate0 + rate1*h)/denominator;
  transition_matrix[1][0] = 1.0 - transition_matrix[1][1];
}


string
CTMarkovModel::tostring() const {
  std::ostringstream oss;
  oss.precision(3);
  oss << "rates  =[" << rate0 << '\t' << rate1 << ']' << endl
      << "U      =[" << U[0][0] << '\t' << U[0][1] << ']' << endl
      << "        [" << U[1][0] << '\t' << U[1][1] << ']' << endl
      << "Uinv   =[" << Uinv[0][0] << '\t' << Uinv[0][1] << ']' << endl
      << "        [" << Uinv[1][0] << '\t' << Uinv[1][1] << ']' << endl
      << "lambda =[" << eigen_values[0] << '\t' << eigen_values[1] << ']';
  return oss.str();
}

std::ostream &
operator<<(std::ostream &os, const CTMarkovModel &ctmm) {
  return os << ctmm.tostring();
}


////////////////////////////////////////////////////////////////////////////////
double
TwoStatesCTMarkovModel::get_trans_prob(const double time_interval,
                                       const size_t start_state,
                                       const size_t end_state) const {
    
  assert(rate0 > 0 && rate1 > 0 && time_interval > 0);
  
  const double h = 1.0 / std::exp(time_interval * (rate0 + rate1));
  const double denominator = rate0 + rate1;
  double prob = 0.0;
  
  if (start_state)
    prob = (rate0 + rate1*h)/denominator;
  else
    prob = (rate0*h + rate1)/denominator;
  
  if (end_state != start_state)
    prob = 1.0 - prob;

  return prob;
}
