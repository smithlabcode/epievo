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

using std::vector;


Environment::Environment(const Path &pa, const Path &pb) {
    bool sa = pa.init_state;
    bool sb = pb.init_state;
    size_t i = 0;
    size_t j = 0;
    tot_time = pa.tot_time;
    while (i < pa.jumps.size() || j < pb.jumps.size()) {
        left.push_back(sa);
        right.push_back(sb);
        if (i < pa.jumps.size() && j < pb.jumps.size()) {
            if (pa.jumps[i] < pb.jumps[j]) {
                breaks.push_back(pa.jumps[i]);
                ++i;
                sa = !sa;
            }
            else if (pa.jumps[i] > pb.jumps[j]) {
                breaks.push_back(pb.jumps[j]);
                ++j;
                sb = !sb;
            }
            else {
                breaks.push_back(pb.jumps[j]);
                ++j;
                ++i;
                sa = !sa;
                sb = !sb;
            }
        }
        else if (i < pa.jumps.size()) {
            breaks.push_back(pa.jumps[i]);
            ++i;
            sa = !sa;
        }
        else {
            breaks.push_back(pb.jumps[j]);
            ++j;
            sb = !sb;
        }
    }
    
    if (breaks.size() == 0 || breaks.back() < tot_time) {
        left.push_back(sa);
        right.push_back(sb);
        breaks.push_back(tot_time);
    }
}
