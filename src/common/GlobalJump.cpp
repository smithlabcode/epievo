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

#include "GlobalJump.hpp"

#include "StateSeq.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <cassert>
#include <sstream>

using std::vector;
using std::string;

static const string ROOT_TAG = "ROOT";
static const string NODE_TAG = "NODE";
static const size_t TAG_LENGTH = 4;

static bool
is_node_line(const string &buffer) {
  return buffer.length() > TAG_LENGTH &&
    (buffer.compare(0, TAG_LENGTH, NODE_TAG) == 0);
}

static bool
is_root_line(const string &buffer) {
  return buffer.length() > TAG_LENGTH &&
    (buffer.compare(0, TAG_LENGTH, ROOT_TAG) == 0);
}

static string
get_node_name(const string &buffer) {
  return string(buffer.begin() + buffer.find(':') + 1, buffer.end());
}

std::istream &
operator>>(std::istream &is, GlobalJump &s) {
  return is >> s.timepoint >> s.position;
}

#include <iomanip>
#include <limits>

std::ostream &
operator<<(std::ostream &os, const GlobalJump &s) {
  return os << std::setprecision(std::numeric_limits<double>::max_digits10)
            << s.timepoint << '\t' << s.position;
}

void
write_root_to_pathfile_global(const string &pathfile, const string &root_name,
                              const StateSeq &root) {
  std::ofstream outpath(pathfile.c_str());
  outpath << ROOT_TAG << ':' << root_name << '\n';
  copy(root.seq.begin(), root.seq.end(),
       std::ostream_iterator<bool>(outpath));
  outpath << '\n';
}

void
append_to_pathfile_global(const string &pathfile, const string &node_name,
                          const vector<GlobalJump> &the_path) {
  std::ofstream outpath(pathfile.c_str(), std::ofstream::app);
  outpath << NODE_TAG << ':' << node_name << '\n';
  for (size_t i = 0; i < the_path.size(); ++i)
    outpath << the_path[i] << '\n';
}

void
read_pathfile_global(const string &pathfile, StateSeq &root,
                     vector<string> &node_names,
                     vector<vector<GlobalJump> > &the_paths) {

  the_paths.clear();

  std::ifstream in(pathfile.c_str());
  if (!in)
    throw std::runtime_error("cannot read: " + pathfile);

  string buffer;
  getline(in, buffer);

  if (!is_root_line(buffer))
    throw std::runtime_error("cannot read root seq: " + pathfile);

  node_names.push_back(get_node_name(buffer));

  getline(in, buffer);
  for (size_t i = 0; i < buffer.length(); ++i)
    root.seq.push_back(buffer[i] == '1');

  the_paths.push_back(vector<GlobalJump>()); // empty on purpose

  std::istringstream iss;
  while (getline(in, buffer)) {
    if (is_node_line(buffer)) {
      node_names.push_back(get_node_name(buffer));
      the_paths.push_back(vector<GlobalJump>());
    }
    else {
      assert(!the_paths.empty());
      iss.clear();
      iss.str(std::move(buffer));
      the_paths.back().push_back(GlobalJump());
      if (!(iss >> the_paths.back().back()))
        throw std::runtime_error("bad line: " + buffer);
    }
  }
}
