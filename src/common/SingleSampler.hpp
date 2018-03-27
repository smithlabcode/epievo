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

void
end_cond_sample(const std::vector<std::vector<double> > &Q,
                const size_t a, const size_t b, const double T,
                std::mt19937 &gen, std::vector<double> &jump_times);



/* Pruning
  - post-order traversal of nodes
  ---- at node v:
  ---- (1) get break points and context, create transition rate matrices
  ---- (2) from bottom up for each break point:
  ------- compute q_k(v), p_j(v)
  ------- (i) if v is a leaf: observed data
  ------- (ii) if v is an internal branching node: two children
  ------- (ii) otherwise: single child


*/
