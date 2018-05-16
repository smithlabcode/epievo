/* Copyright (C) 2018 University of Southern California
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
  TreeHelper(const PhyloTreePreorder &t) : the_tree(t) {
    the_tree.assign_missing_node_names();
    the_tree.get_subtree_sizes(subtree_sizes);
    the_tree.get_node_names(node_names);
    get_parent_id(subtree_sizes, parent_ids);
    the_tree.get_branch_lengths(branches);
    n_nodes = subtree_sizes.size();
  }
  TreeHelper(const double &evo_time) {
    // the_tree is left untouched
    subtree_sizes = {2ul, 1ul};
    node_names = {"root", "leaf"};
    parent_ids = {0, 0};
    branches = {0.0, evo_time};
    n_nodes = 2;
  }

  PhyloTreePreorder the_tree; // tree topology and branch lengths

  /* information we need quick access to, but implicit in the tree */
  std::vector<size_t> subtree_sizes;
  std::vector<std::string> node_names;
  std::vector<size_t> parent_ids;
  std::vector<double> branches;
  size_t n_nodes;
};


#endif
