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

#include <string>
#include <vector>
#include <cassert>
#include <exception>

#include "Segment.hpp"
#include "StateSeq.hpp"

using std::vector;

/* collect rates and interval lengths */
void
collect_segment_info(const vector<double> &rates,
                     const Path &l, const Path &r,
                     vector<SegmentInfo> &seg_info) {
  seg_info.clear();
  seg_info.reserve(l.jumps.size() + r.jumps.size() + 1);
  size_t trip0 = triple2idx(l.init_state, false, r.init_state);
  size_t trip1 = triple2idx(l.init_state, true, r.init_state);
  double prev_time = 0.0;
  vector<double>::const_iterator i(begin(l.jumps)), i_lim = end(l.jumps);
  vector<double>::const_iterator j(begin(r.jumps)), j_lim = end(r.jumps);
  while (i != i_lim && j != j_lim)
    if (*i < *j) {
      seg_info.push_back(SegmentInfo(rates[trip0], rates[trip1], *i-prev_time));
      trip0 = flip_left_bit(trip0);
      trip1 = flip_left_bit(trip1);
      prev_time = *i++;
    }
    else {
      seg_info.push_back(SegmentInfo(rates[trip0], rates[trip1], *j-prev_time));
      trip0 = flip_right_bit(trip0);
      trip1 = flip_right_bit(trip1);
      prev_time = *j++;
    }
  for (; i != i_lim; ++i) {
    seg_info.push_back(SegmentInfo(rates[trip0], rates[trip1], *i - prev_time));
    trip0 = flip_left_bit(trip0);
    trip1 = flip_left_bit(trip1);
    prev_time = *i;
  }
  for (; j != j_lim; ++j) {
    seg_info.push_back(SegmentInfo(rates[trip0], rates[trip1], *j - prev_time));
    trip0 = flip_right_bit(trip0);
    trip1 = flip_right_bit(trip1);
    prev_time = *j;
  }

  assert(l.tot_time == r.tot_time);
  seg_info.push_back(SegmentInfo(rates[trip0], rates[trip1],
                                 l.tot_time - prev_time));
}
