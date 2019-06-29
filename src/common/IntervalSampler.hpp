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

bool
Metropolis_Hastings_interval(const EpiEvoModel &the_model, const TreeHelper &th,
                             const size_t site_id,
                             std::vector<std::vector<Path> > &paths_all_sites_and_branches,
                             std::mt19937 &gen);


/* CASE 1: Path ending in non-bifurcating nodes
  - end-conditioned uniformization sampling
 
   CASE 2: Path involving bifurcating nodes
  - calculate posterior distribution at bifurcating node
  - end-conditioned uniformization sampling along involving intervals

  Compute acceptance rate
*/
