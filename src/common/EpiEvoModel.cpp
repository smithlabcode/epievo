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

#include "EpiEvoModel.hpp"
#include "StateSeq.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <algorithm>

using std::vector;
using std::string;
using std::transform;


string
format_two_by_two(const two_by_two &m) {
  std::ostringstream oss;
  oss << "["
      << std::setw(10) << std::right << m[0][0]
      << std::setw(10) << std::right << m[0][1]
      << "]\n["
      << std::setw(10) << std::right << m[1][0]
      << std::setw(10) << std::right << m[1][1]
      << "]";
  return oss.str();
}

template <class T> string
format_vector(const vector<T> &v) {
  std::ostringstream oss;
  oss << '[';
  if (!v.empty()) {
    oss << v[0];
    for (size_t i = 1; i < v.size(); ++i)
      oss << ',' << v[i];
  }
  oss << ']';
  return oss.str();
}


string
EpiEvoModel::tostring() const {
  std::ostringstream oss;

  oss << "[TREE]" << '\n'
      << t << '\n'
      << "[INIT HORIZ LOG POTENTIALS]\n"
      << format_two_by_two(init_logfac) << '\n'
      << "[INIT HORIZ TRANSITION PROBS]\n"
      << format_two_by_two(init_T) << '\n'
      << "[STATIONARY HORIZ LOG POTENTIALS]\n"
      << format_two_by_two(stationary_logfac) << '\n'
      << "[STATIONARY LOG BASELINE VALUES]\n"
      << format_two_by_two(stationary_logbaseline) << '\n'
      << "[STATIONARY HORIZ TRANSITION PROBS]\n"
      << format_two_by_two(T) << '\n'
      << "[STATIONARY POTENTIALS]\n"
      << format_two_by_two(Q) << '\n'
      << "[TRIPLE RATES]" << '\n';
  for (size_t i = 0; i < triplet_rates.size(); ++i) {
    const size_t left_bit = get_left_bit(i);
    const size_t mid_bit = get_mid_bit(i);
    const size_t right_bit = get_right_bit(i);
    oss << left_bit << mid_bit << right_bit << '\t'
        << triplet_rates[i] << '\n';
  }
  return oss.str();
}


std::ostream &
operator<<(std::ostream &os, const EpiEvoModel &m) {
  return os << m.tostring();
}


void
read_model(const string &param_file, EpiEvoModel &m) {
  std::ifstream in(param_file.c_str());
  if (!in)
    throw std::runtime_error("Could not open file" + param_file);

  string dummy_label;
  in >> dummy_label >> m.t;

  in >> dummy_label;
  assert(dummy_label == "stationary");
  m.stationary_logfac.resize(2, vector<double>(2, 0.0));
  in >> m.stationary_logfac[0][0]
     >> m.stationary_logfac[0][1]
     >> m.stationary_logfac[1][1];
  m.stationary_logfac[1][0] = m.stationary_logfac[0][1];

  in >> dummy_label;
  assert(dummy_label == "baseline");
  m.stationary_logbaseline.resize(2, vector<double>(2, 0.0));
  in >> m.stationary_logbaseline[0][0]
     >> m.stationary_logbaseline[0][1]
     >> m.stationary_logbaseline[1][1];
  m.stationary_logbaseline[1][0] = m.stationary_logbaseline[0][1];

  in >> dummy_label;
  assert(dummy_label == "init");
  m.init_logfac.resize(2, vector<double>(2, 0.0));
  in >> m.init_logfac[0][0]
     >> m.init_logfac[0][1]
     >> m.init_logfac[1][1];
  m.init_logfac[1][0] = m.init_logfac[0][1];

  m.initialize();
}


void
EpiEvoModel::get_rates(vector<double> &rates) const {
  rates.resize(8, 0.0);
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      for (size_t k = 0; k < 2; ++k)
        rates[triple2idx(i, j, k)] = std::exp(stationary_logfac[i][1-j] +
                                              stationary_logfac[1-j][k] +
                                              stationary_logbaseline[i][k]);
}


////////////////////////////////////////////////////////////////////////////////
// stationary and Markov chain
////////////////////////////////////////////////////////////////////////////////

void
potential_to_transition_prob(const two_by_two &Q, two_by_two &T) {

  const double delta = sqrt(pow(Q[0][0] - Q[1][1], 2) + 4*Q[0][1]*Q[1][0]);

  // transition probability matrix
  T = Q;

  // compute the diagonal entries
  const double diag_denom = Q[0][0] + Q[1][1] + delta;
  T[1][1] = 2*Q[1][1]/diag_denom;
  T[0][0] = 2*Q[0][0]/diag_denom;

  // now compute the anti-diagonal entries
  const double anti_numer = 4*Q[0][1]*Q[1][0];
  T[0][1] = anti_numer/(pow(Q[0][0] + delta, 2) - Q[1][1]*Q[1][1]);
  T[1][0] = anti_numer/(pow(Q[1][1] + delta, 2) - Q[0][0]*Q[0][0]);
}


void
scale_rates(const vector<double> &rates, const vector<double> &branches,
            vector<double> &scaled_rates, vector<double> &scaled_branches) {

  // pair-wise potentials Q
  two_by_two rate_mat(2, vector<double>(2, 1.0));
  rate_mat[0][0] = rate_mat[0][1] * sqrt(rates[2]/rates[0]);
  rate_mat[1][1] = rate_mat[0][1] * sqrt(rates[2]/rates[0]) * (rates[1]/rates[3]);

  // transition probability matrix
  two_by_two trans_mat;
  potential_to_transition_prob(rate_mat, trans_mat);

  vector<double> pi(2, 0);
  pi[1] = (1.0 - trans_mat[0][0])/(2.0 - trans_mat[0][0] - trans_mat[1][1]);
  pi[0] = (1.0 - pi[1]);
  vector<double> stationary_prob(rates.size(), 0.0);

  double unit = 0.0;
  for (size_t i = 0; i < rates.size(); ++i) {
    const size_t left = get_left_bit(i);
    const size_t mid = get_mid_bit(i);
    const size_t right = get_right_bit(i);

    stationary_prob[i] = pi[left]*trans_mat[left][mid]*trans_mat[mid][right];

    unit += stationary_prob[i]*rates[i];
  }

  scaled_rates = rates;
  transform(scaled_rates.begin(), scaled_rates.end(),
            scaled_rates.begin(), std::bind2nd(std::divides<double>(), unit));

  scaled_branches = branches;
  transform(scaled_branches.begin(), scaled_branches.end(),
            scaled_branches.begin(),
            std::bind2nd(std::multiplies<double>(), unit));
}


/* This function takes the parameter values for the model, which were
   most likely read from a parameter file, and computes the initial
   horizontal transitions (init_T), the stationary horizontal
   transitions (T), the stationary rates (Q), and the rates
   for triples */
void
EpiEvoModel::initialize() {

  // convert the initial potentials into transition probs
  two_by_two init_Q(2, vector<double>(2, 0.0));
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      init_Q[i][j] = exp(init_logfac[i][j]);
  potential_to_transition_prob(init_Q, init_T);

  // convert the stationary potentials into transition probs
  Q = two_by_two(2, vector<double>(2, 0.0));
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      Q[i][j] = exp(stationary_logfac[i][j]);
  potential_to_transition_prob(Q, T);

  // get the 1D vector of rates for the triples
  triplet_rates.resize(8, 0.0);
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      for (size_t k = 0; k < 2; ++k)
        triplet_rates[triple2idx(i, j, k)] =
          std::exp(stationary_logfac[i][1-j] +
                   stationary_logfac[1-j][k] +
                   stationary_logbaseline[i][k]);
}


void
sample_state_sequence(const size_t n_sites, const two_by_two &trans_prob,
                      std::mt19937 &gen, vector<char> &sequence) {

  sequence = vector<char>(n_sites, true);

  const double pi1 =
    (1.0 - trans_prob[0][0])/(2.0 - trans_prob[1][1] - trans_prob[0][0]);

  std::uniform_real_distribution<double> unif(0.0, 1.0);
  sequence[0] = (unif(gen) < pi1);

  for (size_t i = 1; i < n_sites; ++i) {
    const double r = unif(gen);
    const double p = sequence[i - 1] ? trans_prob[1][1] : trans_prob[0][0];
    sequence[i] = ((r <= p) ? sequence[i - 1] : !sequence[i - 1]);
  }
}


void
EpiEvoModel::sample_state_sequence_init(const size_t n_sites, std::mt19937 &gen,
                                        vector<char> &sequence) const {
  sample_state_sequence(n_sites, init_T, gen, sequence);
}

void
EpiEvoModel::sample_state_sequence_stationary(const size_t n_sites,
                                              std::mt19937 &gen,
                                              vector<char> &sequence) const {
  sample_state_sequence(n_sites, T, gen, sequence);
}
