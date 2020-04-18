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

#include <vector>

#include "Path.hpp"
#include "TreeHelper.hpp"
#include "epievo_utils.hpp"

struct FelsHelper {
  FelsHelper() {}
  FelsHelper(const std::vector<double> &_p,
             const std::vector<double> &_q) : p(_p), q(_q) {}
  std::vector<double> p;
  std::vector<double> q;
};

void
expectation_sufficient_statistics(const std::vector<double> &rates,
                                  const TreeHelper &th,
                                  const std::vector<std::vector<Path> > &paths,
                                  std::vector<std::vector<double> > &J,
                                  std::vector<std::vector<double> > &D);


void
update_paths_indep(const std::vector<double> &rates, const TreeHelper &th,
                   std::vector<std::vector<Path> > &paths,
                   std::mt19937 &gen);


void
compute_sufficient_statistics(const std::vector<std::vector<Path> > &paths,
                              std::vector<std::vector<double> > &J,
                              std::vector<std::vector<double> > &D);


void
estimate_rates_indep(const std::vector<std::vector<double> > &J,
                     const std::vector<std::vector<double> > &D,
                     std::vector<double> &rates,
                     TreeHelper &th);


void
estimate_rates_and_branches_indep(const std::vector<std::vector<double> > &J,
                                  const std::vector<std::vector<double> > &D,
                                  std::vector<double> &rates,
                                  TreeHelper &th,
                                  std::vector<std::vector<Path> > &paths);

void
estimate_root_stationary(const std::vector<std::vector<Path> > &paths,
                         std::vector<double> &pi);
