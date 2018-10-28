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

#include <string>
#include <vector>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <random>
#include <functional>
#include <iomanip>

#include "SingleSiteSampler.hpp"
#include "PhyloTreePreorder.hpp"
#include "Path.hpp"

#include "StateSeq.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::pair;
using std::make_pair;

using std::bind;
using std::placeholders::_1;
using std::multiplies;
using std::runtime_error;

// collect rates and interval lengths
static void
collect_segment_info(const vector<double> &triplet_rates,
                     const Path &l, const Path &r,
                     vector<SegmentInfo> &seg_info, size_t &n_segs) {
  Environment env(l, r);
  const size_t n_intervals = env.left.size();
  seg_info = vector<SegmentInfo>(n_intervals);

  for (size_t i = 0; i < n_intervals; ++i) {
    const size_t pattern0 = triple2idx(env.left[i], false, env.right[i]);
    const size_t pattern1 = triple2idx(env.left[i], true, env.right[i]);
    seg_info[i] = SegmentInfo(triplet_rates[pattern0], triplet_rates[pattern1],
                              env.breaks[i] - (i == 0 ? 0.0 : env.breaks[i-1]));
    assert(seg_info[i].len > 0.0);
    n_segs++;
  }
}


bool
Metropolis_Hastings_site(const EpiEvoModel &the_model, const TreeHelper &th,
                         const size_t site_id, vector<vector<Path> > &paths,
                         std::mt19937 &gen, vector<Path> &proposed_path,
                         const bool always_accept) {
  // get rates and lengths each interval [seg_info: node x site]
  // and label all the nodes (including induced breakpoints)
  size_t n_segs = 0;
  vector<vector<SegmentInfo> > seg_info(th.n_nodes);
  vector<size_t> all_node_ids;
  for (size_t node_id = 1; node_id < th.n_nodes; ++node_id) {
    collect_segment_info(the_model.triplet_rates,
                         paths[node_id][site_id - 1],
                         paths[node_id][site_id + 1], seg_info[node_id],
                         n_segs);
  }

  std::uniform_int_distribution<size_t> seg_sampler(0, n_segs - 1);
  const size_t seg_to_update = seg_sampler(gen);
  
  // locate the target segment
  size_t seg_scanner = 0;
  for (size_t node_id = 1; node_id < th.n_nodes; ++node_id) {
    size_t n_intervals = seg_info[node_id - 1];
    for (size_t i = 0; i < n_intervals; ++i) {
      seg_scanner++;
      if (seg_scanner == seg_to_update) {
        // the target segment located
        // need to know whether it's joining other branches
        if (i > 0 && i < n_intervals - 1)
          propose_h_interval(th, site_id, paths, seg_info[node_id][i]);
        else if (i == 0)
          if (th.is_root(th.parent_ids[node_id]))
            propose_root_interval(th, site_id, paths, seg_info[node_id][i]);
          else
            propose_branching_intervals(th, site_id, paths, seg_info, node_id,
                                        i);
        else if (i == n_intervals - 1)
          if (th.is_leaf(node_id))
            propose_root_interval(th, site_id, paths, seg_info[node_id][i]);
          else
            propose_branching_intervals(th, site_id, paths, seg_info, node_id,
                                        i);
      }

  bool accepted = false;
    
  return accepted;
}
