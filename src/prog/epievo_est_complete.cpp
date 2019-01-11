/* Copyright (C) 2018 University of Southern California
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

    string outfile;
    bool VERBOSE = false;
    bool OPTBRANCH = false;
    bool scale_the_rates = true;
    string param_file;
    string tree_file;

    ////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]), "estimate parameters from"
                           " complete data (site-specific paths)",
                           "<path-file>");
    opt_parse.add_opt("param", 'p', "initial parameter file",
                      true, param_file);
    opt_parse.add_opt("tree", 't', "initial tree file in newick format",
                      true, tree_file);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("no-rate-scaling", 'S', "do not scale rates",
                      false, scale_the_rates);
    opt_parse.add_opt("branch", 'b', "optimize branch lengths as well",
                      false, OPTBRANCH);
    opt_parse.add_opt("output", 'o', "output parameter file",
                      false, outfile);
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string path_file(leftover_args.front());
    ////////////////////////////////////////////////////////////////////////

    if (VERBOSE)
      cerr << "[READING PATHS FILE: " << path_file << "]" << endl;
    vector<vector<Path> > all_paths; // along multiple branches
    vector<string> node_names;
    read_paths(path_file, node_names, all_paths);
    const size_t n_sites = all_paths.back().size();

    if (VERBOSE)
      cerr << "[READING PARAMETER FILE: " << param_file << endl;
    EpiEvoModel the_model;
    read_model(param_file, the_model);
    if (VERBOSE)
      cerr << the_model << endl;

    if (VERBOSE)
      cerr << "[READING TREE FILE: " << tree_file << "]" << endl;
    PhyloTreePreorder the_tree; // tree topology and branch lengths
    std::ifstream tree_in(tree_file.c_str());
    if (!tree_in || !(tree_in >> the_tree))
      throw std::runtime_error("bad tree file: " + tree_file);
    const size_t n_nodes = the_tree.get_size();
    TreeHelper th(the_tree);

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
    else
      compute_estimates_rates_and_branches(VERBOSE, param_tol, all_paths,
                                           th, the_model);

    estimate_root_distribution(all_paths, the_model);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    if (!out)
      throw std::runtime_error("bad output file: " + outfile);

    out << the_model.format_for_param_file() << endl;
    the_tree.set_branch_lengths(th.branches);
    out << the_tree << endl;
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
