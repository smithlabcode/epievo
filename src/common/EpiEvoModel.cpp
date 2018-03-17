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
#include <bitset>

using std::vector;
using std::string;
using std::transform;
using std::bitset;


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
      << "[STATIONARY HORIZ TRANSITION PROBS]\n"
      << format_two_by_two(T) << '\n'
      << "[STATIONARY LOG BASELINE VALUES]\n"
      << format_two_by_two(stationary_logbaseline) << '\n'
      << "[STATIONARY POTENTIALS]\n"
      << format_two_by_two(Q) << '\n'
      << "[TRIPLE RATES]" << '\n';
  oss << bitset<3>(0) << '\t' << triplet_rates[0];
  for (size_t i = 1; i < triplet_rates.size(); ++i)
    oss << '\n' << bitset<3>(i) << '\t' << triplet_rates[i];
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

  /* read the phylogenetic tree */
  string dummy_label;
  in >> dummy_label >> m.t;

  /* read the stationary distribution */
  in >> dummy_label;
  assert(dummy_label == "stationary");
  m.stationary_logfac.resize(2, vector<double>(2, 0.0));
  in >> m.stationary_logfac[0][0]
     >> m.stationary_logfac[0][1]
     >> m.stationary_logfac[1][1];
  m.stationary_logfac[1][0] = m.stationary_logfac[0][1];

  /* read the baseline */
  in >> dummy_label;
  assert(dummy_label == "baseline");
  m.stationary_logbaseline.resize(2, vector<double>(2, 0.0));
  in >> m.stationary_logbaseline[0][0]
     >> m.stationary_logbaseline[0][1]
     >> m.stationary_logbaseline[1][1];
  m.stationary_logbaseline[1][0] = m.stationary_logbaseline[0][1];

  /* read the initial distribution (at root) */
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
compute_stationary_state_proportions(const two_by_two &T, vector<double> &pi) {
  pi.resize(2);
  pi[1] = (1.0 - T[0][0])/(2.0 - T[0][0] - T[1][1]);
  pi[0] = (1.0 - pi[1]);
}

void
EpiEvoModel::get_stationary_state_proportions(vector<double> &pi) const {
  compute_stationary_state_proportions(T, pi);
}

void
compute_stationary_triplet_proportions(const two_by_two &T,
                                       vector<double> &props) {
  vector<double> pi;
  compute_stationary_state_proportions(T, pi);
  props.resize(EpiEvoModel::n_triplets);
  for (size_t i = 0; i < EpiEvoModel::n_triplets; ++i) {
    const size_t left = get_left_bit(i);
    const size_t mid = get_mid_bit(i);
    const size_t right = get_right_bit(i);
    props[i] = pi[left]*T[left][mid]*T[mid][right]; // T is member variable
  }
}


void
EpiEvoModel::get_stationary_triplet_proportions(vector<double> &props) const {
  compute_stationary_triplet_proportions(T, props);
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

  // proportions of triplets at stationarity
  vector<double> stationary_triplet_props;
  compute_stationary_triplet_proportions(trans_mat, stationary_triplet_props);

  double unit = 0.0;
  for (size_t i = 0; i < stationary_triplet_props.size(); ++i)
    unit += stationary_triplet_props[i]*rates[i];

  scaled_rates = rates;
  transform(scaled_rates.begin(), scaled_rates.end(),
            scaled_rates.begin(), std::bind2nd(std::divides<double>(), unit));

  scaled_branches = branches;
  transform(scaled_branches.begin(), scaled_branches.end(),
            scaled_branches.begin(),
            std::bind2nd(std::multiplies<double>(), unit));
}

void
EpiEvoModel::compute_triplet_rates() {
  triplet_rates.resize(n_triplets, 0.0);
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      for (size_t k = 0; k < 2; ++k)
        triplet_rates[triple2idx(i, j, k)] =
          std::exp(stationary_logfac[i][1-j] +
                   stationary_logfac[1-j][k] +
                   stationary_logbaseline[i][k]);
}


/* This function takes the parameter values for the model, which were
   most likely read from a parameter file, and computes the initial
   horizontal transitions (init_T), the stationary horizontal
   transitions (T), the stationary rates (Q), and the rates
   for triples */
void
EpiEvoModel::initialize() {

  // make sure every node has a name ADS: modify below so that all
  // nodes always have names; check for this rather than attempt to
  // fix it
  t.assign_missing_node_names();
  t.get_subtree_sizes(subtree_sizes);
  t.get_node_names(node_names);
  get_parent_id(subtree_sizes, parent_ids);
  t.get_branch_lengths(branches);

  // convert the initial potentials into transition probs
  two_by_two init_Q(2, vector<double>(2, 0.0));
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      init_Q[i][j] = exp(init_logfac[i][j]);
  potential_to_transition_prob(init_Q, init_T);

  // initialize stationary potentials
  Q = two_by_two(2, vector<double>(2, 0.0));
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      Q[i][j] = exp(stationary_logfac[i][j]);

  // convert the stationary potentials into transition probs
  potential_to_transition_prob(Q, T);

  // get the 1D vector of rates for the triples
  compute_triplet_rates();
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

double
EpiEvoModel::substitutions_per_site(const vector<double> &triplet_props) const {
  assert(triplet_props.size() == n_triplets);
  double total_rate = 0;
  for (size_t i = 0; i < n_triplets; ++i)
    total_rate += triplet_rates[i]*triplet_props[i];
  return total_rate;
}


void
pairwise_potentials_from_triplet_rates(const vector<double> &triplet_rates,
                                       two_by_two &Q) {
  // pair-wise potentials Q from rates lambda_ijk
  Q = vector<vector<double> >(2, vector<double>(2, 1.0));
  const double death_birth_ratio = triplet_rates[2]/triplet_rates[0];
  const double expand_contract_ratio = triplet_rates[1]/triplet_rates[3];
  Q[0][0] = Q[0][1] * sqrt(death_birth_ratio);
  Q[1][1] = Q[0][1] * sqrt(death_birth_ratio) * expand_contract_ratio;
}


void
stationary_from_horiz_trans_prob(const two_by_two &T, vector<double> &pi) {
  pi = vector<double>(2, 0.0);
  pi[1] = (1.0 - T[0][0])/(2.0 - T[0][0] - T[1][1]);
  pi[0] = 1.0 - pi[1];
}

double
rate_scaling_factor(const vector<double> &pi,
                    const two_by_two &T,
                    const vector<double> &triplet_rates) {

  double mu_rate_value = 0.0;
  for (size_t i = 0; i < triplet_rates.size(); ++i) {
    const size_t l = get_left_bit(i);
    const size_t m = get_mid_bit(i);
    const size_t r = get_right_bit(i);
    mu_rate_value += (pi[l]*T[l][m]*T[m][r]) * triplet_rates[i];
  }
  return mu_rate_value;
}

double
rate_scaling_factor(const vector<double> &triplet_rates) {

  // pair-wise potentials Q from lambda_ijk
  two_by_two Q;
  pairwise_potentials_from_triplet_rates(triplet_rates, Q);

  // transition probability matrix T from Q
  two_by_two T;
  potential_to_transition_prob(Q, T);

  // stationary rates pi from T
  vector<double> pi;
  stationary_from_horiz_trans_prob(T, pi);

  return rate_scaling_factor(pi, T, triplet_rates);
}
