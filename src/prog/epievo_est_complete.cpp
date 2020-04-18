/* Copyright (C) 2019 University of Southern California
 *                    Jianghan Qu, Andrew D Smith and Xiaojing Ji
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
#include <iostream>
#include <fstream>
#include <cmath>   /* exp, sqrt, pow fabs*/
#include <numeric>  /* std::inner_product, accumulate*/
#include <functional>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "ParamEstimation.hpp"
#include "PhyloTreePreorder.hpp"
#include "Path.hpp"
#include "EpiEvoModel.hpp"
#include "TreeHelper.hpp"

#include "epievo_utils.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::ofstream;
using std::runtime_error;

int main(int argc, const char **argv) {

  try {

    static const double param_tol = 1e-10;

    bool VERBOSE = false;
    bool single_branch = false;
    bool optimize_branches = false;

    string outfile;
    string tree_file, treefile_updated;

    ////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]), "estimate parameters from"
                           " complete data (site-specific paths)",
                           "<param> (<treefile>) <path_file>");
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("branch", 'b', "optimize branch lengths as well",
                      false, optimize_branches);
    opt_parse.add_opt("single_branch", 'T', "pairwise process (assumes no tree)",
                      false, single_branch);
    opt_parse.add_opt("output", 'o', "output parameter file",
                      true, outfile);
    opt_parse.add_opt("outtree", 't', "output file of tree",
                      false, treefile_updated);
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
    if (leftover_args.size() == 2) {
      if (!single_branch) {
        cerr << opt_parse.help_message() << endl;
        return EXIT_SUCCESS;
      }
    }
    else if (leftover_args.size() != 3) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    else {
      tree_file = leftover_args[1];
    }
    const string param_file(leftover_args.front());
    const string path_file(leftover_args.back());
    ////////////////////////////////////////////////////////////////////////

    /* LOADING PATHS */
    if (VERBOSE)
      cerr << "[READING PATHS FILE: " << path_file << "]" << endl;
    vector<string> node_names;
    vector<vector<Path> > orig_paths; // [nodes] x [sites]
    vector<vector<Path> > paths; // [sites] x [nodes]
    read_paths(path_file, node_names, orig_paths);

    // transposing
    const size_t n_sites = orig_paths[1].size();
    const size_t n_nodes = orig_paths.size();
    paths.resize(n_sites);
    for (size_t i = 0; i < n_sites; ++i) {
      paths[i].resize(n_nodes);
      for (size_t b = 1; b < n_nodes; ++b)
        paths[i][b] = orig_paths[b][i];
    }

    /* LOADING (FAKE) TREE */
    PhyloTreePreorder the_tree; // tree topology and branch lengths
    TreeHelper th;
    if (single_branch) {
      if (VERBOSE)
        cerr << "[INITIALIZING TWO NODE TREE WITH TIME: "
        << paths.front().back().tot_time << "]" << endl;
      th = TreeHelper(paths.front().back().tot_time);
    }
    else {
      if (VERBOSE)
        cerr << "[READING TREE: " << tree_file << "]" << endl;
      std::ifstream tree_in(tree_file);
      if (!tree_in || !(tree_in >> the_tree))
        throw runtime_error("bad tree file: " + tree_file);
      th = TreeHelper(the_tree);
    }

    if (VERBOSE)
      cerr << "[READING PARAMETER FILE: " << param_file << endl;
    EpiEvoModel the_model;
    read_model(param_file, the_model);
    the_model.scale_triplet_rates();
    if (VERBOSE)
      cerr << the_model << endl;

    if (VERBOSE) {
      cerr << "[INITIAL GUESS AT BRANCH LENGTHS]" << endl;
      for (size_t i = 0; i < th.node_names.size(); ++i)
        cerr << th.node_names[i] << '\t' << th.branches[i] << endl;
    }

    if (VERBOSE)
      cerr << "n_nodes=" << n_nodes << endl
           << "n_sites=" << n_sites << endl;

    if (!optimize_branches) {
      estimate_rates(VERBOSE, param_tol, paths, the_model);
    }
    else {
      estimate_rates_and_branches(VERBOSE, param_tol, paths, th, the_model);
      scale_jump_times(paths, th);
      the_tree.set_branch_lengths(th.branches);
    }

    if (VERBOSE)
      cerr << "[WRITING PARAMETERS]\n" <<  the_model << endl;
    ofstream out(outfile);
    if (!out)
      throw runtime_error("bad output file: " + outfile);
    out << the_model.format_for_param_file() << endl;

    if (optimize_branches && !treefile_updated.empty()) {
      ofstream out_tree(treefile_updated);
      if (!out_tree)
        throw runtime_error("bad output param file: " + treefile_updated);
      out_tree << the_tree << endl;
    }

    if (VERBOSE)
      cerr << "[FINISHED]" << endl;
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
