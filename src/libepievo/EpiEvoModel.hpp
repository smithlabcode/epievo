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

#ifndef EPIEVOMODEL_HPP
#define EPIEVOMODEL_HPP

#include <string>
#include <vector>
#include <random>
#include <array>

#include "epievo_utils.hpp"

struct EpiEvoModel {

  // constant for the number of triplets (avoiding magic)
  static const size_t n_triplets = 8;

  two_by_two stationary_baseline; // symmetric part in J-P model
  two_by_two T;       // horizontal transition probs (stationary)
  two_by_two Q;       // pair-wise potential densities (stationary)

  std::array<double, n_triplets> triplet_rates;
  //std::vector<double> triplet_rates; // rates for triples

  // mutators
  void rebuild_from_triplet_rates(const std::array<double, n_triplets> &updated_rates);
  void scale_triplet_rates();
  void initialize(); // assuming rates are always scaled to one change
                     // per site per unit time

  // accessors
  double substitutions_per_site(const std::vector<double> &triplet_props) const;
  bool is_unit_rate() const;
  void get_stationary_state_proportions(std::vector<double> &pi) const;
  void get_stationary_triplet_proportions(std::vector<double> &props) const;

  void sample_state_sequence(const size_t n_sites, std::mt19937 &gen,
                             std::vector<bool> &sequence) const;

  std::string tostring() const;
  std::string format_for_param_file() const;

private:
  void compute_triplet_rates();
};


std::ostream &
operator<<(std::ostream &os, const EpiEvoModel &m);

void
read_model(const bool SCALE, const std::string &param_file, EpiEvoModel &m);

void
read_model(const std::string &param_file, EpiEvoModel &m);

double
rate_scaling_factor(const std::array<double, 8> &triplet_rates);


void
decompose(const std::vector<double> &rates,
          std::vector<double> &eigen_vals,
          std::vector<std::vector<double> > &U,
          std::vector<std::vector<double> > &Uinv);

#endif
