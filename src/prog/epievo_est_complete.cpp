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
#include "StateSeq.hpp"
#include "EpiEvoModel.hpp"
#include "TreeHelper.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;


int main(int argc, const char **argv) {

  try {

    static const double param_tol = 1e-10;

    bool VERBOSE = false;
    bool OPTBRANCH = false;
    bool ONEBRANCH = false;

    string outfile;
    string treefile_updated;

    ////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]), "estimate parameters from"
                           " complete data (site-specific paths)",
                           "<param> <treefile> <path_file>");
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("branch", 'b', "optimize branch lengths as well",
                      false, OPTBRANCH);
    opt_parse.add_opt("one-branch", 'T', "one-branch tree", false,
                      ONEBRANCH);
    opt_parse.add_opt("output", 'o',
                      "output parameter file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("outtree", 't',
                      "output file of tree (default: stdout)",
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
    if (leftover_args.size() < 3) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string param_file(leftover_args[0]);
    const string tree_file(leftover_args[1]);
    const string path_file(leftover_args[2]);
    ////////////////////////////////////////////////////////////////////////

    if (VERBOSE)
      cerr << "[READING PATHS FILE: " << path_file << "]" << endl;
    vector<vector<Path> > all_paths; // along multiple branches
    vector<string> node_names;
    read_paths(path_file, node_names, all_paths);
    const size_t n_sites = all_paths.back().size();

    /* LOADING (FAKE) TREE */
    if (VERBOSE)
      cerr << "[READING TREE: " << tree_file << "]" << endl;
    PhyloTreePreorder the_tree; // tree topology and branch lengths
    TreeHelper th;
    if (ONEBRANCH) {
      if (VERBOSE)
        cerr << "initializing two node tree with time: " << tree_file << endl;
      th = TreeHelper(std::stod(tree_file));
    } else {
      cerr << "reading tree file: " << tree_file << endl;
      std::ifstream tree_in(tree_file.c_str());
      if (!tree_in || !(tree_in >> the_tree))
        throw std::runtime_error("bad tree file: " + tree_file);
      th = TreeHelper(the_tree);
    }
    const size_t n_nodes = the_tree.get_size();

    if (VERBOSE)
      cerr << "[READING PARAMETER FILE: " << param_file << endl;
    EpiEvoModel the_model;
    read_model(param_file, the_model);
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

    if (!OPTBRANCH) {
      compute_estimates_for_rates_only(VERBOSE, param_tol,
                                       all_paths, the_model);
    }
    else {
      compute_estimates_rates_and_branches(VERBOSE, param_tol, all_paths,
                                           th, the_model);
      scale_jump_times(all_paths, th);
      the_tree.set_branch_lengths(th.branches);
    }

    estimate_root_distribution(all_paths, the_model);

    if (VERBOSE)
      cerr << "[WRITING PARAMETERS]\n" <<  the_model << endl;
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    if (!out)
      throw std::runtime_error("bad output file: " + outfile);

    out << the_model.format_for_param_file() << endl;

    if (OPTBRANCH) {
      std::ofstream of_tree;
      if (!treefile_updated.empty()) of_tree.open(treefile_updated.c_str());
      std::ostream out_tree(treefile_updated.empty() ?
                            std::cout.rdbuf() : of_tree.rdbuf());
      if (!out_tree)
        throw std::runtime_error("bad output param file: " + treefile_updated);
      out_tree << the_tree << endl;
    }

    if (VERBOSE)
      cerr << "[FINISHED.]" << endl;
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
