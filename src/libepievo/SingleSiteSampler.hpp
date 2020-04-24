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

#ifndef SINGLE_SITE_SAMPLER_HPP
#define SINGLE_SITE_SAMPLER_HPP

#include <vector>

#include "Path.hpp"
#include "Segment.hpp"
#include "EpiEvoModel.hpp"
#include "TreeHelper.hpp"

struct FelsHelper;

class SingleSiteSampler {
public:

  SingleSiteSampler(const size_t n_burn_in, const size_t n_batch);

  void
  reset(const EpiEvoModel &the_model, std::vector<std::vector<Path> > &paths);

  bool
  Metropolis_Hastings_site(const EpiEvoModel &the_model, const TreeHelper &th,
                           const size_t site_id,
                           std::vector<std::vector<Path> > &paths,
                           std::mt19937 &gen);

  void
  run_mcmc(const EpiEvoModel &the_model, const TreeHelper &th,
           std::vector<std::vector<Path> > &paths,
           std::mt19937 &gen,
           std::vector<std::vector<double> > &J_all_sites,
           std::vector<std::vector<double> > &D_all_sites,
           double &acceptance_rate);

  // MCMC parameter constants
  bool SAMPLE_ROOT;
  size_t burn_in;
  size_t batch;

private:

  size_t
  single_iteration(const EpiEvoModel &the_model, const TreeHelper &th,
                   std::vector<std::vector<Path> > &paths,
                   std::mt19937 &gen);

  // internal variables for faster computation
  std::vector<double> tri_llh;
  double log_rates[8];
  std::uniform_real_distribution<double> unif;

  // scratch variables
  std::vector<double> J;
  std::vector<double> D;
  std::vector<Path> proposed_path;
};

#endif
