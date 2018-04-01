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

#include <vector>

#include "Path.hpp"
#include "EpiEvoModel.hpp" /* model_param */

void
end_cond_sample(const double rate0, const double rate1,
                const size_t a, const size_t b, const double T,
                std::mt19937 &gen, std::vector<double> &jump_times);


void
pruning(const std::vector<double> &triplet_rates,
        const std::vector<size_t> &subtree_sizes,
        const size_t site,
        const std::vector<std::vector<Path> > &all_paths,
        const std::vector<std::vector<std::vector<double> > > &all_interval_rates,
        const std::vector<std::vector<double> > & all_interval_lengths,
        std::vector<std::vector<std::vector<double> > > &all_p);

void
downward_sampling(const std::vector<double> &triplet_rates,
                  const std::vector<size_t> &subtree_sizes,
                  const size_t site,
                  const std::vector<std::vector<Path> > &all_paths,
                  const std::vector<std::vector<double> > &root_trans_prob,
                  const std::vector<std::vector<std::vector<double> > > &all_p,
                  std::mt19937 &gen,
                  std::vector<Path> &new_path);


void
gibbs_site(const EpiEvoModel &the_model,
           const size_t site,
           std::vector<std::vector<Path> > &all_paths,
           std::mt19937 &gen,
           std::vector<Path> &new_path);

/* Pruning
  - post-order traversal of nodes
  ---- at node v:
  ---- (1) get break points and context, create transition rate matrices
  ---- (2) from bottom up for each break point:
  ------- compute q_k(v), p_j(v)
  ------- (i) if v is a leaf: observed data
  ------- (ii) if v is an internal branching node: two children
  ------- (ii) otherwise: single child

  Downward sampling
  - preorder traversal
  - at root, compute posterior probability of 0 vs 1
  - at each internal node and break points within branch do conditional posterior sampling

  Compute acceptance rate (todo)
*/
