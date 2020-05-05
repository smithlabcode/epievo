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

#ifndef SEGMENT_HPP
#define SEGMENT_HPP

#include <string>
#include <vector>
#include <array>
#include <cassert>

#include "Path.hpp"

struct SegmentInfo {
  SegmentInfo() {}
  SegmentInfo(const double r0, const double r1,
              const size_t i0, const size_t i1, const double l) :
    rate0(r0), rate1(r1), trip0(i0), trip1(i1), len(l) {}
  double rate0;
  double rate1;
  size_t trip0;
  size_t trip1;
  double len;
};

void
collect_segment_info(const std::array<double, 8> &rates,
                     const Path &l, const Path &r,
                     std::vector<SegmentInfo> &seg_info);
#endif
