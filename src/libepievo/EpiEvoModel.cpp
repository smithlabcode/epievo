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

#include "EpiEvoModel.hpp"
#include "epievo_utils.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <functional>
#include <bitset>

using std::vector;
using std::string;
using std::transform;
using std::bitset;
using std::bind;
using std::placeholders::_1;

using std::cerr;
using std::endl;

static bool
is_probability_distribution(const vector<double> &v) {
  return std::all_of(begin(v), end(v), bind(std::greater<double>(), _1, 0.0)) &&
    std::accumulate(begin(v), end(v), 0.0) == 1.0;
}


string
format_two_by_two(const two_by_two &m) {
  std::ostringstream oss;
  oss << "["
      << std::setw(10) << std::right << m(0, 0)
      << std::setw(10) << std::right << m(0, 1)
      << "]\n["
      << std::setw(10) << std::right << m(1, 0)
      << std::setw(10) << std::right << m(1, 1)
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


void
horiz_trans_prob_to_horiz_potential(const two_by_two &T,
                                    two_by_two &Q) {

  assert(is_probability_distribution(T[0]) &&
         is_probability_distribution(T[1]));

  // Gibbs measure pair-wise potentials
  Q = T;
  Q(0, 0) = 1.0 - T(0, 1);
  Q(0, 1) = sqrt(T(0, 1)*T(1, 0));
  Q(1, 0) = Q(0, 1);
  Q(1, 1) = 1.0 - T(1, 0);
}

/* This function assumes phi(0,1) = 0, and will produce Q that can be
 * used for some computations, but will not always recover the Q that
 * might be expected.
 */
void
triplet_rates_to_horiz_potential(const vector<double> &triplet_rates,
                                 two_by_two &Q) {
  // pair-wise potentials Q from rates lambda_ijk
  Q = two_by_two (1.0, 1.0, 1.0, 1.0);
  const double death_birth_ratio = triplet_rates[2]/triplet_rates[0];
  const double expand_contract_ratio = triplet_rates[1]/triplet_rates[3];
  Q(0, 0) = Q(0, 1) * sqrt(death_birth_ratio);
  Q(1, 1) = Q(0, 1) * sqrt(death_birth_ratio) * expand_contract_ratio;
}


void
horiz_potential_to_horiz_trans_prob(const two_by_two &Q,
                                    two_by_two &T) {

  const double delta = sqrt(pow(Q(0, 0) - Q(1, 1), 2) + 4*Q(0, 1)*Q(1, 0));

  // transition probability matrix
  T = Q;

  // compute the diagonal entries
  const double diag_denom = Q(0, 0) + Q(1, 1) + delta;
  T(1, 1) = 2*Q(1, 1)/diag_denom;
  T(0, 0) = 2*Q(0, 0)/diag_denom;

  // now compute the anti-diagonal entries
  T(0, 1) = 1.0 - T(0, 0);
  T(1, 0) = 1.0 - T(1, 1);

  assert(is_probability_distribution(T[0]) &&
         is_probability_distribution(T[1]));
}


void
triplet_rates_to_horiz_trans_prob(const vector<double> &triplet_rates,
                                  two_by_two &T) {

  // first compute approximation to pair-wise potentials Q from the triplet rates
  two_by_two Q_proportional;
  triplet_rates_to_horiz_potential(triplet_rates, Q_proportional);

  // now use that approx. to get the exact T
  horiz_potential_to_horiz_trans_prob(Q_proportional, T);
}



void
horiz_trans_prob_to_horiz_stationary(const two_by_two &T,
                                     vector<double> &pi) {
  pi = vector<double>(2, 0.0);
  pi[1] = (1.0 - T(0, 0))/(2.0 - T(0, 0) - T(1, 1));
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
    mu_rate_value += pi[l]*T(l, m)*T(m, r)*triplet_rates[i];
  }
  return mu_rate_value;
}


double
rate_scaling_factor(const vector<double> &triplet_rates) {

  // first compute approximation to pair-wise potentials Q from the triplet rates
  two_by_two Q_proportional;
  triplet_rates_to_horiz_potential(triplet_rates, Q_proportional);

  // now use that approx. to get the exact T
  two_by_two T;
  horiz_potential_to_horiz_trans_prob(Q_proportional, T);

  // stationary rates pi from T
  vector<double> pi;
  horiz_trans_prob_to_horiz_stationary(T, pi);

  return rate_scaling_factor(pi, T, triplet_rates);
}


std::string
EpiEvoModel::format_for_param_file() const {
  std::ostringstream oss;
  oss << "stationary\t" << T(0, 0) << '\t' << T(1, 1) << endl
      << "baseline\t"
      << stationary_baseline(0, 0) << '\t'
      << stationary_baseline(1, 1) << endl
      << "init\t" << init_T(0, 0) << '\t' << init_T(1, 1);
  return oss.str();
}


string
EpiEvoModel::tostring() const {
  std::ostringstream oss;

  oss << "[INIT HORIZ TRANSITION PROBS]\n"
      << format_two_by_two(init_T) << '\n'
      << "[STATIONARY HORIZ TRANSITION PROBS]\n"
      << format_two_by_two(T) << '\n'
      << "[STATIONARY BASELINE VALUES]\n"
      << format_two_by_two(stationary_baseline) << '\n'
      << "[STATIONARY POTENTIALS DENSITIES]\n"
      << format_two_by_two(Q) << '\n'
      << "[TRIPLE RATES]\n";
  oss << bitset<3>(0) << '\t' << triplet_rates[0];
  for (size_t i = 1; i < triplet_rates.size(); ++i)
    oss << '\n' << bitset<3>(i) << '\t' << triplet_rates[i];

  const double mu = rate_scaling_factor(triplet_rates);
  oss << "\n[UNIT TIME TRANSITIONS: " << mu << "]";

  return oss.str();
}


std::ostream &
operator<<(std::ostream &os, const EpiEvoModel &m) {
  return os << m.tostring();
}


void
compute_stationary_state_proportions(const two_by_two &T,
                                     vector<double> &pi) {
  pi.resize(2);
  pi[1] = (1.0 - T(0, 0))/(2.0 - T(0, 0) - T(1, 1));
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
    props[i] = pi[left]*T(left, mid)*T(mid, right); // T is member variable
  }
}


void
EpiEvoModel::get_stationary_triplet_proportions(vector<double> &props) const {
  compute_stationary_triplet_proportions(T, props);
}


void
scale_rates(const vector<double> &rates, const vector<double> &branches,
            vector<double> &scaled_rates, vector<double> &scaled_branches) {

  double unit = rate_scaling_factor(rates);

  scaled_rates = rates;
  transform(scaled_rates.begin(), scaled_rates.end(), scaled_rates.begin(),
            std::bind(std::divides<double>(), std::placeholders::_1, unit));

  scaled_branches = branches;
  transform(scaled_branches.begin(), scaled_branches.end(),
            scaled_branches.begin(),
            std::bind(std::multiplies<double>(), std::placeholders::_1, unit));
}


void
sample_state_sequence(const size_t n_sites,
                      const two_by_two &trans_prob,
                      std::mt19937 &gen, vector<bool> &sequence) {

  sequence = vector<bool>(n_sites, true);

  const double pi1 =
    (1.0 - trans_prob(0, 0))/(2.0 - trans_prob(1, 1) - trans_prob(0, 0));

  std::uniform_real_distribution<double> unif(0.0, 1.0);
  sequence[0] = (unif(gen) < pi1);
  for (size_t i = 1; i < n_sites; ++i) {
    const double r = unif(gen);
    const double p = sequence[i - 1] ? trans_prob(1, 1) : trans_prob(0, 0);
    sequence[i] = ((r <= p) ? sequence[i - 1] : (!sequence[i - 1]));
  }
}


void
EpiEvoModel::sample_state_sequence_init(const size_t n_sites, std::mt19937 &gen,
                                        vector<bool> &sequence) const {
  sample_state_sequence(n_sites, init_T, gen, sequence);
}

void
EpiEvoModel::sample_state_sequence_stationary(const size_t n_sites,
                                              std::mt19937 &gen,
                                              vector<bool> &sequence) const {
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


bool
EpiEvoModel::is_unit_rate() const {
  return rate_scaling_factor(triplet_rates) == 1.0;
}


/* read the model file
 */
void
read_model(const string &param_file, EpiEvoModel &m) {

  std::ifstream in(param_file.c_str());
  if (!in)
    throw std::runtime_error("Could not open file: " + param_file);

  string dummy_label;
  in >> dummy_label;
  if(dummy_label == "stationary") {
    /* read the stationary distribution */
    m.T.reset();
    in >> m.T(0, 0) >> m.T(1, 1);
    m.T(1, 0) = 1.0 - m.T(1, 1);
    m.T(0, 1) = 1.0 - m.T(0, 0);

    assert(is_probability_distribution(m.T[0]) &&
           is_probability_distribution(m.T[1]));

    /* read the baseline */
    in >> dummy_label;
    assert(dummy_label == "baseline");
    m.stationary_baseline.reset();
    in >> m.stationary_baseline(0, 0)
    >> m.stationary_baseline(1, 1);

    /* read the initial distribution (at root) */
    in >> dummy_label;
    assert(dummy_label == "init");
    m.init_T.reset();
    in >> m.init_T(0, 0) >> m.init_T(1, 1);
    m.init_T(1, 0) = 1.0 - m.init_T(1, 1);
    m.init_T(0, 1) = 1.0 - m.init_T(0, 0);

    assert(is_probability_distribution(m.init_T[0]) &&
           is_probability_distribution(m.init_T[1]));

    m.initialize();
  }
  else {
    assert(dummy_label == "000");
    m.T.reset();
    m.stationary_baseline.reset();
    m.init_T.reset();

    /* read triplet transition rates */
    vector<double> rates;
    rates.resize(8, 0.0);
    in >> rates[0];
    for (size_t i = 1; i < 8; i++)
      in >> dummy_label >> rates[i];

    /* take care of constraints between rates */
    // lambda_100 = lambda_001
    rates[4] = rates[1];
    // lambda_110 = lambda_011
    rates[6] = rates[3];
    // B * C^2 * M = D * E^2 * S
    // load B,D,E,C,M, compute S
    rates[7] = (rates[0] * rates[6] * rates[6] * rates[5]) /
    (rates[2] * rates[4] * rates[4]);

    m.rebuild_from_triplet_rates(rates);
    m.init_T = m.T;
  }
}

void
EpiEvoModel::scale_triplet_rates() {
  const double mu = rate_scaling_factor(triplet_rates);
  for (size_t i = 0; i < n_triplets; ++i)
    triplet_rates[i] /= mu;
}

/* This function takes the parameter values for the model, which were
   most likely read from a parameter file, and computes the initial
   horizontal transitions (init_T), the stationary horizontal
   transitions (T), the stationary rates (Q), and the rates
   for triples */
void
EpiEvoModel::initialize() {

  // convert target transition probs into Gibbs pair-wise potentials
  Q = two_by_two();
  horiz_trans_prob_to_horiz_potential(T, Q);

  // get the 1D vector of rates for the triples
  compute_triplet_rates();

  // scale for one change per site per unit time
  // scale_triplet_rates();
}


void
EpiEvoModel::compute_triplet_rates() {
  triplet_rates.resize(n_triplets, 0.0);

  double exp_stat_bl = exp(stationary_baseline(0,0));
  triplet_rates[triple2idx(0,0,0)] = Q(0,1)*Q(1,0)*exp_stat_bl;
  triplet_rates[triple2idx(0,1,0)] = Q(0,0)*Q(0,0)*exp_stat_bl;

  exp_stat_bl = exp(stationary_baseline(0,1));
  triplet_rates[triple2idx(0,0,1)] = Q(0,1)*Q(1,1)*exp_stat_bl;
  triplet_rates[triple2idx(0,1,1)] = Q(0,0)*Q(0,1)*exp_stat_bl;

  exp_stat_bl = exp(stationary_baseline(1, 0));
  triplet_rates[triple2idx(1,0,0)] = Q(1,1)*Q(1,0)*exp_stat_bl;
  triplet_rates[triple2idx(1,1,0)] = Q(1,0)*Q(0,0)*exp_stat_bl;

  exp_stat_bl = exp(stationary_baseline(1, 1));
  triplet_rates[triple2idx(1,0,1)] = Q(1,1)*Q(1,1)*exp_stat_bl;
  triplet_rates[triple2idx(1,1,1)] = Q(1,0)*Q(0,1)*exp_stat_bl;

  // for (size_t i = 0; i < 2; ++i)
  //   for (size_t j = 0; j < 2; ++j)
  //     for (size_t k = 0; k < 2; ++k)
  //       triplet_rates[triple2idx(i, j, k)] =
  //         Q(i, 1-j)*Q(1-j, k)*std::exp(stationary_baseline(i, k));
}


void
EpiEvoModel::rebuild_from_triplet_rates(const vector<double> &updated_rates) {
  assert(updated_rates.size() == n_triplets);
  assert(updated_rates[1] == updated_rates[4]);
  assert(updated_rates[3] == updated_rates[6]);

  triplet_rates = updated_rates;

  // recompute T using the updated triplet rates
  triplet_rates_to_horiz_trans_prob(triplet_rates, T);

  // recompute Q using the updated T
  horiz_trans_prob_to_horiz_potential(T, Q);

  // recompute stationary log baseline
  const double log_Q01 = log(Q(0, 1));
  const double log_Q10 = log(Q(1, 0));
  const double log_Q11 = log(Q(1, 1));
  stationary_baseline(0, 0) = log(triplet_rates[0]) - (log_Q01 + log_Q10);
  stationary_baseline(0, 1) = log(triplet_rates[1]) - (log_Q01 + log_Q11);
  stationary_baseline(1, 0) = log(triplet_rates[4]) - (log_Q11 + log_Q10);
  stationary_baseline(1, 1) = log(triplet_rates[7]) - (log_Q10 + log_Q01);

  const double centering_point = stationary_baseline(0, 1);
  stationary_baseline(0, 0) -= centering_point;
  stationary_baseline(0, 1) -= centering_point;
  stationary_baseline(1, 0) -= centering_point;
  stationary_baseline(1, 1) -= centering_point;

  assert(stationary_baseline(0, 1) == stationary_baseline(1, 0));
}


void
decompose(const vector<double> &rates, // rate0 and rate1
          vector<double> &eigen_vals,
          two_by_two &U, two_by_two &Uinv) {

  // Q = U*D*Uinv
  const double sum_rate = rates[0] + rates[1];
  eigen_vals = {0.0, -sum_rate};

  U(0, 0) = 1.0;
  U(0, 1) = rates[0];
  U(1, 0) = 1.0;
  U(1, 1) = -rates[1];

  Uinv(0, 0) = rates[1]/sum_rate;
  Uinv(0, 1) = rates[0]/sum_rate;
  Uinv(1, 0) = 1.0/sum_rate;
  Uinv(1, 1) = -1.0/sum_rate;
}
