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
#include <iostream>
#include <fstream>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "PhyloTreePreorder.hpp"

#include "GlobalJump.hpp"
#include "StateSeq.hpp"
#include "path.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::runtime_error;


static void
read_states_file(const string &statesfile,
                 vector<string> &node_names,
                 vector<StateSeq> &the_states) {

  std::ifstream in(statesfile.c_str());
  if (!in)
    throw runtime_error("cannot read states file: " + statesfile);

  string buffer;
  getline(in, buffer);
  if (buffer[0] == '#') // ADS: probably should adopt some convention
    buffer = buffer.substr(1);

  std::istringstream iss(buffer);
  string tmp_name;
  while (iss >> tmp_name)
    node_names.push_back(tmp_name);

  const size_t n_nodes = node_names.size();

  size_t the_site = 0; // dummy
  the_states.resize(n_nodes);
  while (getline(in, buffer)) {
    iss.clear();
    iss.str(std::move(buffer));

    iss >> the_site;

    size_t val_count = 0;
    char tmp_value = 0;
    while (val_count < n_nodes && iss >> tmp_value)
      the_states[val_count++].seq.push_back(tmp_value == '1');

    if (val_count != n_nodes)
      throw runtime_error("bad line in states file");
  }
}


static void
assign_changes_to_sites(const vector<GlobalJump> &global_path,
                        vector<Path> &by_site) {

  const size_t n_changes = global_path.size();
  for (size_t i = 0; i < n_changes; ++i)
    by_site[global_path[i].position].jumps.push_back(global_path[i].timepoint);
}


int main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;

    OptionParser opt_parse(strip_path(argv[0]), "convert path file format",
                           "<treefile> <statefile> <pathfile> <outfile>");
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 4) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string treefile(leftover_args[0]);
    const string statesfile(leftover_args[1]);
    const string pathsfile(leftover_args[2]);
    const string outfile(leftover_args[3]);
    ////////////////////////////////////////////////////////////////////////

    if (VERBOSE)
      cerr << "[reading paths: " << pathsfile << "]" << endl;

    PhyloTreePreorder the_tree; // tree topology and branch lengths
    std::ifstream tree_in(treefile.c_str());
    if (!tree_in || !(tree_in >> the_tree))
      throw runtime_error("cannot read tree file: " + treefile);

    vector<size_t> subtree_sizes;
    the_tree.get_subtree_sizes(subtree_sizes);
    vector<string> node_names;
    the_tree.get_node_names(node_names);
    vector<size_t> parent_ids;
    get_parent_id(subtree_sizes, parent_ids);
    vector<double> branch_lengths;
    the_tree.get_branch_lengths(branch_lengths);

    const size_t n_nodes = node_names.size();

    StateSeq root;
    vector<string> node_names_from_pathsfile;
    vector<vector<GlobalJump> > the_paths; // along multiple branches
    read_pathfile_global(pathsfile, root, node_names_from_pathsfile, the_paths);

    write_root_to_pathfile_global(outfile, node_names.front(), root);

    vector<StateSeq> the_states;
    vector<string> node_names_from_statesfile;
    read_states_file(statesfile, node_names_from_statesfile, the_states);

    assert(node_names == node_names_from_statesfile &&
           node_names == node_names_from_pathsfile);

    const size_t n_sites = root.seq.size();

    std::ofstream out(outfile.c_str());
    if (!out)
      throw std::runtime_error("bad output file: " + outfile);

    for (size_t node_id = 1; node_id < n_nodes; ++node_id) {

      const size_t parent_id = parent_ids[node_id];
      const double branch_len = branch_lengths[node_id];

      vector<Path> path_by_site(n_sites);
      for (size_t j = 0; j < path_by_site.size(); ++j) {
        path_by_site[j].init_state = the_states[parent_id].seq[j];
        path_by_site[j].tot_time = branch_len;
      }
      assign_changes_to_sites(the_paths[node_id], path_by_site);

      out << "NODE:" << node_names[node_id] << endl;
      for (size_t i = 0; i < path_by_site.size(); ++i)
        out << i << '\t' << path_by_site[i] << endl;
    }
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
