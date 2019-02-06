/* Copyright (C) 2019 University of Southern California
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

#include <iostream>
#include <string>
#include <vector>

class StateSeq;

struct GlobalJump {
  /* GlobalJump: represents a change in the boolean state within a
     sequence of boolean states. The identity of the change is deduced
     from the identity prior to the change, which ultimately is
     deduced from some ground state (i.e. the root sequence).
   */
  GlobalJump() {}
  GlobalJump(const double tp, const size_t pos) : timepoint(tp), position(pos) {}
  double timepoint; // time of change
  size_t position; // position of change
  bool operator<(const GlobalJump &other) const {
    return timepoint < other.timepoint;
  }
};

std::ostream &
operator<<(std::ostream &os, const GlobalJump &s);

std::istream &
operator>>(std::istream &is, GlobalJump &s);

void
write_root_to_pathfile_global(const std::string &pathfile,
                              const std::string &root_name,
                              const StateSeq &root);

void
append_to_pathfile_global(const std::string &pathfile,
                          const std::string &node_name,
                          const std::vector<GlobalJump> &the_path);

void
read_pathfile_global(const std::string &pathfile, StateSeq &root,
                     std::vector<std::string> &node_names,
                     std::vector<std::vector<GlobalJump> > &the_paths);
