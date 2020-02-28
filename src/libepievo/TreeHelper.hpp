/* Copyright (C) 2019 University of Southern California
 *                    Jianghan Qu, Xiaojing Ji and Andrew D Smith
 *
 * Author: Andrew D. Smith
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

/* This class exists to keep together the auxiliary information used
   to work with the tree, and to ensure everything is consistent. This
   functionality could be included in PhyloTreePreorder, because that
   child class of PhyloTree exists to support the functionality we
   need here, but I'm doing it this way for convenience. This class
   should be kept small...
 */

#ifndef TREEHELPER_HPP
#define TREEHELPER_HPP

#include <string>
#include <vector>

#include "PhyloTreePreorder.hpp"

struct TreeHelper {

  TreeHelper() {}
  TreeHelper(const PhyloTreePreorder &t);
  TreeHelper(const double &evo_time);

  PhyloTreePreorder the_tree; // tree topology and branch lengths

  /* information we need quick access to, but implicit in the tree */
  std::vector<size_t> subtree_sizes;
  std::vector<std::string> node_names;
  std::vector<size_t> parent_ids;
  std::vector<double> branches;
  size_t n_nodes;

  // FUNCTIONS BELOW HERE
  bool is_leaf(const size_t node_id) const;
  bool is_root(const size_t node_id) const;

};

class ChildSet {
public:
  ChildSet(const std::vector<size_t> &v, const size_t ni) :
    offset(v.begin() + ni), node_id(ni), ch_id(1) {}

  size_t operator*() const {return node_id + ch_id;}

  bool good() const {return ch_id < *offset;}

  ChildSet & operator++() {
    ch_id += *(offset + ch_id);
    return *this;
  }
  ChildSet operator++(int) { // this one is slower...
    ChildSet tmp(*this);
    operator++();
    return tmp;
  }

private:
  const std::vector<size_t>::const_iterator offset;
  size_t node_id;
  size_t ch_id;
};

#endif
