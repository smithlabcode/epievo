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
#include "Path.hpp"
#include "TripletSampler.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::runtime_error;


static void
process_branch_above_node(const double report_interval,
                          const string &outfile,
                          const StateSeq &initial_states,
                          const vector<GlobalJump> &the_jumps) {

  std::ofstream out(outfile.c_str());
  if (!out)
    throw runtime_error("cannot open output file: " + outfile);

  TripletSampler ts(initial_states);

  double current_time = 0.0;
  double next_report_time = 0.0;

  const size_t n_changes = the_jumps.size();
  for (size_t i = 0; i < n_changes; ++i) {
    const double next_jump_time = the_jumps[i].timepoint;

    if (next_report_time < next_jump_time) {
      StateSeq s;
      ts.get_sequence(s);
      // in case there are no changes in some reporting interval
      while (next_report_time < next_jump_time) {
        cout << next_report_time << '\t' << current_time << '\t'
             << next_jump_time << '\t'
             << current_time + next_jump_time << endl;
        out << s << endl;
        next_report_time += report_interval;
      }
    }
    ts.mutate(the_jumps[i].position);
    current_time = next_jump_time;
  }
}


int main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    double report_interval = 0.0;
    size_t n_reports = 0;

    OptionParser opt_parse(strip_path(argv[0]), "convert path file format",
                           "<node-name> <treefile> <statefile> "
                           "<pathfile> <outfile>");
    opt_parse.add_opt("interval", 'i', "time interval for extracting states",
                      false, report_interval);
    opt_parse.add_opt("reports", 'r', "number of time-points to report",
                      false, n_reports);
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
    if ((n_reports == 0 && report_interval == 0.0) ||
        (n_reports > 0  && report_interval >  0.0)) {
      cerr << "exactly one of options \'i\' and \'r\' required" << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 5) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string node_name(leftover_args[0]);
    const string treefile(leftover_args[1]);
    const string statesfile(leftover_args[2]);
    const string pathsfile(leftover_args[3]);
    const string outfile(leftover_args[4]);
    ////////////////////////////////////////////////////////////////////////

    if (VERBOSE)
      cerr << "[reading tree: " << treefile << "]" << endl;
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

    if (VERBOSE)
      cerr << "[reading paths: " << pathsfile << "]" << endl;
    StateSeq root;
    vector<string> node_names_from_pathsfile;
    vector<vector<GlobalJump> > the_paths; // along multiple branches
    read_pathfile_global(pathsfile, root, node_names_from_pathsfile, the_paths);

    if (VERBOSE)
      cerr << "[reading states: " << statesfile << "]" << endl;
    vector<StateSeq> the_states;
    vector<string> node_names_from_statesfile;
    read_states_file(statesfile, node_names_from_statesfile, the_states);

    assert(node_names == node_names_from_statesfile &&
           node_names == node_names_from_pathsfile);

    // get the id for the desired node
    vector<string>::const_iterator name_idx =
      find(begin(node_names), end(node_names), node_name);
    if (name_idx == node_names.end())
      throw runtime_error("invalid node name: " + node_name);

    const size_t node_id = name_idx - node_names.begin();
    const size_t parent_id = parent_ids[node_id];
    const double branch_length = branch_lengths[node_id];

    if (n_reports > 0)
      report_interval = branch_length/n_reports;

    if (VERBOSE) {
      cerr << "node name: " << node_name << endl
           << "node id: " << node_id << endl
           << "parent name: " << node_names[parent_id] << endl
           << "parent id: " << parent_id << endl
           << "branch length: " << branch_length << endl
           << "report interval: " << report_interval << endl
           << "number of reports: " << n_reports << endl
           << "total jumps: " << the_paths[node_id].size() << endl;
    }

    // add a fake jump to the end of the jump sequence above the
    // desired node; this will never coincide with a reporting time,
    // so the effect will not be observed.
    the_paths[node_id].push_back(GlobalJump(branch_lengths[node_id], 1));

    process_branch_above_node(report_interval, outfile,
                              the_states[parent_id], the_paths[node_id]);
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
