/* Copyright (C) 2018 University of Southern California
 *                    Jianghan Qu and Andrew D Smith
 *
 * Author: Jianghan Qu and Andrew D. Smith
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

#include "smithlab_utils.hpp"

using std::vector;
using std::string;

struct Path {

  Path() : init_state(false), tot_time(0.0) {}

  bool init_state;
  double tot_time;
  vector<double> jumps;

  bool state_at_time(const double t) const;
  bool end_state() const {
    return (jumps.size() % 2 == 0) ? init_state : !init_state;
  }
};

std::ostream &
operator<<(std::ostream &os, const Path &p);

void to_path(const bool s, const string jumps, Path &p);


void initialize_paths(const std::vector<bool> &seq, const double tot_time,
                      std::vector<Path> &paths);

void read_paths(const string path_file, vector<vector<Path> > &paths);

void
read_paths(const string &path_file,
           std::vector<std::string> &node_names,
           std::vector<std::vector<Path> > &paths);

void get_initial_seq(const vector<Path> &paths, vector<bool> &seq);

void end_sequence(const vector<Path> &paths, vector<bool> &seq);

void sequence_at_time(const vector<Path> &paths, const double t,
                      vector<bool> &seq);

////////////////////////////////////////////////////////////////////////////////

struct Environment {
  // states on the left
  vector<bool> left;
  // states on the right
  vector<bool> right;
  // time points of environment state changes, including tot_time
  vector<double> breaks;
  double tot_time;

  Environment(const Path &pa, const Path &pb);
};



////////////////////////////////////////////////////////////////////////////////
struct TriplePath {
  vector<size_t> states; /*triplet states, length k*/
  vector<double> breaks; /*start with first jump, end with total_time, length k*/
  vector<size_t> jump_context_freq; /*context frequency of jumps at middle site*/

  TriplePath(const Path &l, const Path &m, const Path &r);
  void time_by_context(vector<double> &tbc) const;
};

////////////////////////////////////////////////////////////////////////////////
struct PathContextStat {
  vector<double> jumps_in_context;
  vector<double> time_in_context;

  PathContextStat(const Path &l, const Path &m, const Path &r);
};

#endif
