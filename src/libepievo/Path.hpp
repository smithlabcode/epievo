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

#ifndef PATH_HPP
#define PATH_HPP

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cassert>

#include "epievo_utils.hpp"

struct Path {

  Path() : init_state(false), tot_time(0.0) {}
  Path(const bool is, const double tt) :
    init_state(is), tot_time(tt), jumps(std::vector<double>()) {}
  Path(const bool is, const double tt, const std::vector<double> &j) :
    init_state(is), tot_time(tt), jumps(j) {}

  bool init_state;
  double tot_time;
  std::vector<double> jumps;

  bool is_valid() const {
    for (size_t i = 1; i < jumps.size(); ++i)
      if (jumps[i] < jumps[i-1])
        return false;
    return true;
  }

  bool state_at_time(const double t) const;
  bool end_state() const {
    return (jumps.size() % 2 == 0) ? init_state : !init_state;
  }
  void scale_to_unit_length();
};

std::ostream &
operator<<(std::ostream &os, const Path &p);

void
to_path(const bool s, const std::string &jumps, Path &p);

void
initialize_paths(const std::vector<bool> &seq, const double tot_time,
                 std::vector<Path> &paths);

void
read_paths(const std::string &path_file, std::vector<std::vector<Path> > &paths);

void
read_paths(const std::string &path_file,
           std::vector<std::string> &node_names,
           std::vector<std::vector<Path> > &paths);

void
get_seq_init(const std::vector<Path> &paths, state_seq &seq);

void
get_seq_end(const std::vector<Path> &paths, state_seq &seq);

void
get_seq_at_time(const double t, const std::vector<Path> &paths,
                state_seq &seq);

void
add_sufficient_statistics(const Path &left, const Path &mid, const Path &right,
                          std::vector<double> &J, std::vector<double> &D);

#endif
