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

// template <typename T>
// struct two_by_two {
//   two_by_two() : x00(0), x01(0), x10(0), x11(0) {}
//   T x00;
//   T x01;
//   T x10;
//   T x11;

//   template <class U> T
//   operator()(const U first, const U second) const {
//     if (first) {
//       if (second) return x11;
//       else return x10;
//     }
//     else {
//       if (second) return x01;
//       else return x00;
//     }
//   }

//   template <class U> T &
//   operator()(const U first, const U second) {
//     if (first) {
//       if (second) return x11;
//       else return x10;
//     }
//     else {
//       if (second) return x01;
//       else return x00;
//     }
//   }
// };

typedef std::vector<std::vector<double> > two_by_two;

struct EpiEvoModel {

  PhyloTreePreorder t; // tree topology and branch lengths

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

  string tostring() const;

private:
  void get_rates(std::vector<double> &rates) const;
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

#endif
