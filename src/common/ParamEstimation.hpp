/* Copyright (C) 2018 University of Southern California
 *                    Jianghan Qu, Andrew D Smith and Xiaojing Ji
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

#include <vector>

void
compute_estimates_for_rates_only(const bool VERBOSE,
                                 const double param_tol,
                                 const std::vector<double> &J,
                                 const std::vector<double> &D,
                                 EpiEvoModel &the_model);

void
compute_estimates_for_rates_only(const bool VERBOSE,
                                 const double param_tol,
                                 const std::vector<std::vector<Path> > &all_paths,
                                 EpiEvoModel &the_model);

void
compute_estimates_rates_and_branches(const bool VERBOSE,
                                     const double param_tol,
                                     std::vector<std::vector<Path> > &all_paths,
                                     TreeHelper &th,
                                     EpiEvoModel &the_model);

void
estimate_root_distribution(const std::vector<std::vector<Path> > &all_paths,
                           EpiEvoModel &the_model);


void
get_sufficient_statistics(const std::vector<std::vector<Path> > &all_paths,
                          std::vector<double> &J, std::vector<double> &D);

void
estimate_sufficient_statistics_by_simulation(const std::vector<std::vector<Path> > &all_paths,
                                             const EpiEvoModel the_model,
                                             const TreeHelper &th,
                                             std::vector<double> &J,
                                             std::vector<double> &D,
                                             std::mt19937 &gen,
                                             const size_t n = 10000);

    
#endif
