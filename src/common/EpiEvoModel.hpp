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

#include "PhyloTreePreorder.hpp"

typedef std::vector<std::vector<double> > two_by_two;

struct EpiEvoModel {

  PhyloTreePreorder t; // tree topology and branch lengths

  /* information we need quick access to, but implicit in the tree */
  std::vector<size_t> subtree_sizes;
  std::vector<std::string> node_names;
  std::vector<size_t> parent_ids;
  std::vector<double> branches;

  // "logfac" means log of potential factors
  two_by_two stationary_logfac;      // log of stationary potentials for pairs
  two_by_two stationary_logbaseline; // symmetric part in J-P model
  two_by_two init_logfac;            // log of initial potentials for pairs

  two_by_two T;       // horizontal transition probs (stationary)
  two_by_two init_T;  // horizontal transition probs (initial)
  two_by_two Q;       // pair-wise potentials (stationary)

  std::vector<double> triplet_rates; // rates for triples

  void initialize();
  void sample_state_sequence_init(const size_t n_sites, std::mt19937 &gen,
                                  std::vector<char> &sequence) const;
  void sample_state_sequence_stationary(const size_t n_sites, std::mt19937 &gen,
                                        std::vector<char> &sequence) const;

  std::string tostring() const;

  std::string root_name() const {return node_names[0];}

private:
  void compute_triplet_rates();
};

std::ostream &
operator<<(std::ostream &os, const EpiEvoModel &m);

void
read_model(const std::string &param_file, EpiEvoModel &m);

void
potential_to_transition_prob(const two_by_two &Q, two_by_two &T);

void
scale_rates(const std::vector<double> &rates,
            const std::vector<double> &branches,
            std::vector<double> &scaled_rates,
            std::vector<double> &scaled_branches);

double
rate_scaling_factor(const std::vector<double> &triplet_rates);

#endif
