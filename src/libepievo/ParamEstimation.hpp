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

#ifndef PARAMETER_ESTIMATION_HPP
#define PARAMETER_ESTIMATION_HPP

#include "ContinuousTimeMarkovModel.hpp"
#include "Path.hpp"
#include "TreeHelper.hpp"
#include "EpiEvoModel.hpp"

#include "epievo_utils.hpp"

#include <vector>
#include <array>

double
estimate_rates(const bool VERBOSE,
               const double param_tol,
               const std::vector<std::vector<double> > &J,
               const std::vector<std::vector<double> > &D,
               EpiEvoModel &the_model);

double
estimate_rates(const bool VERBOSE,
               const double param_tol,
               const std::vector<std::vector<Path> > &all_paths,
               EpiEvoModel &the_model);


double
estimate_rates_and_branches(const bool VERBOSE,
                            const double param_tol,
                            const std::vector<std::vector<double> > &J,
                            const std::vector<std::vector<double> > &D,
                            TreeHelper &th,
                            EpiEvoModel &the_model);

double
estimate_rates_and_branches(const bool VERBOSE, const double param_tol,
                            std::vector<std::vector<Path> > &all_paths,
                            TreeHelper &th, EpiEvoModel &the_model);

// get sufficient statistics from the whole tree
void
get_sufficient_statistics(const std::vector<std::vector<Path> > &all_paths,
                          std::vector<double> &J, std::vector<double> &D);

// get sufficient statistics from each branch
void
get_sufficient_statistics(const std::vector<std::vector<Path> > &all_paths,
                          std::vector<std::vector<double> > &J,
                          std::vector<std::vector<double> > &D);

// count duplet frequencies from root sequence
void
get_root_frequencies(const std::vector<std::vector<Path> > &all_paths,
                     two_by_two &counts);

// scale jump times according to updated branch lengths
void
scale_jump_times(std::vector<std::vector<Path> > &all_paths,
                 const TreeHelper &th);

void
set_one_change_per_site_per_unit_time(std::array<double, 8> &rates,
                                      std::vector<double> &branches);

#endif
