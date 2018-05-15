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

#ifndef EPIEVOMODEL_HPP
#define EPIEVOMODEL_HPP

#include <string>
#include <vector>
#include <random>

typedef std::vector<std::vector<double> > two_by_two;

struct EpiEvoModel {

  void
  rebuild_from_triplet_rates(const std::vector<double> &triplet_rates);

  two_by_two stationary_logbaseline; // symmetric part in J-P model
  two_by_two T;       // horizontal transition probs (stationary)
  two_by_two init_T;  // horizontal transition probs (initial)
  two_by_two Q;       // pair-wise potentials (stationary)

  std::vector<double> triplet_rates; // rates for triples

  double substitutions_per_site(const std::vector<double> &triplet_props) const;
  void get_stationary_state_proportions(std::vector<double> &pi) const;
  void get_stationary_triplet_proportions(std::vector<double> &props) const;

  void scale_triplet_rates();
  void initialize(); // assuming rates are always scaled to one change
                     // per site per unit time

  void sample_state_sequence_init(const size_t n_sites, std::mt19937 &gen,
                                  std::vector<char> &sequence) const;
  void sample_state_sequence_stationary(const size_t n_sites, std::mt19937 &gen,
                                        std::vector<char> &sequence) const;

  std::string tostring() const;
  std::string format_for_param_file() const;

  static const size_t n_triplets = 8;

private:
  void compute_triplet_rates();
};


std::ostream &
operator<<(std::ostream &os, const EpiEvoModel &m);

void
read_model(const bool SCALE, const std::string &param_file, EpiEvoModel &m);

void
read_model(const std::string &param_file, EpiEvoModel &m);

void
scale_rates(const std::vector<double> &rates,
            const std::vector<double> &branches,
            std::vector<double> &scaled_rates,
            std::vector<double> &scaled_branches);

double
rate_scaling_factor(const std::vector<double> &triplet_rates);

void
continuous_time_trans_prob_mat(const double rate0, const double rate1,
                               const double time_interval,
                               std::vector<std::vector<double> >
                               &transition_matrix);

void
decompose(const std::vector<double> &rates,
          std::vector<double> &eigen_vals,
          std::vector<std::vector<double> > &U,
          std::vector<std::vector<double> > &Uinv);

#endif
