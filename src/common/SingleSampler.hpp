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
#include "EpiEvoModel.hpp"
#include "TreeHelper.hpp"
#include "EndCondSampling.hpp"

struct SegmentInfo {
  SegmentInfo() {}
  SegmentInfo(const double r0, const double r1, const double l) :
    rate0(r0), rate1(r1), len(l) {}
  double rate0;
  double rate1;
  double len;
};

struct FelsHelper {
  FelsHelper() {}
  FelsHelper(const std::vector<std::vector<double> > &_p,
             const std::vector<double> &_q) : p(_p), q(_q) {}
  std::vector<std::vector<double> > p;
  std::vector<double> q;
};

void
pruning(const TreeHelper &th, const size_t site_id,
        const std::vector<std::vector<Path> > &paths_all_sites_and_branches,
        const std::vector<std::vector<SegmentInfo> > &seg_info,
        std::vector<FelsHelper> &fh);

void
downward_sampling(const TreeHelper &th, const size_t site_id,
                  const std::vector<std::vector<Path> > &paths_all_sites_and_branches,
                  const std::vector<std::vector<double> > &root_trans_prob,
                  const std::vector<std::vector<SegmentInfo> > &seg_info,
                  const std::vector<FelsHelper> &fh,
                  std::mt19937 &gen,
                  std::vector<Path> &new_path);

void
Metropolis_Hastings_site(const EpiEvoModel &the_model, const TreeHelper &th,
                         const size_t site_id,
                         std::vector<std::vector<Path> > &paths_all_sites_and_branches,
                         std::mt19937 &gen,
                         std::vector<Path> &sampled_path);

void
Gibbs_site(const EpiEvoModel &the_model, const TreeHelper &th,
           const size_t site_id,
           std::vector<std::vector<Path> > &paths_all_sites_and_branches,
           std::mt19937 &gen,
           std::vector<Path> &sampled_path);


void
add_sufficient_statistics(const Path &left, const Path &mid, const Path &right,
                          std::vector<double> &J, std::vector<double> &D);

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
