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


static void
delete_last_branch(vector<vector<Path> > &paths, vector<string> &node_names,
                   TreeHelper &th) {
  paths.pop_back();
  node_names.pop_back();
  th.branches.pop_back();
  th.parent_ids.pop_back();
  th.subtree_sizes.pop_back();
  th.subtree_sizes[0]--;
  th.n_nodes--;
}


int main(int argc, const char **argv) {

  try {

    static const double param_tol = 1e-10;

    string outfile;
    bool VERBOSE = false;
    bool OPTBRANCH = false;
    bool scale_the_rates = true;
    string param_file;
    string tree_file;
    
    size_t rng_seed = std::numeric_limits<size_t>::max();
    size_t iteration = 10;

    ////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]), "estimate parameters from"
                           " initial data (site-specific paths)",
                           "<path-file>");
    opt_parse.add_opt("param", 'p', "initial parameter file",
                      true, param_file);
    opt_parse.add_opt("tree", 't', "initial tree file in newick format",
                      true, tree_file);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);
    opt_parse.add_opt("no-rate-scaling", 'S', "do not scale rates",
                      false, scale_the_rates);
    opt_parse.add_opt("branch", 'b', "optimize branch lengths as well",
                      false, OPTBRANCH);
    opt_parse.add_opt("output", 'o', "output parameter file",
                      false, outfile);
    opt_parse.add_opt("iteration", 'i', "number of MCMC-EM iterations",
                      false, iteration);
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
    TreeHelper th(the_tree);
    // remove unwanted branch for now
    delete_last_branch(all_paths, node_names, th);
    size_t n_nodes = all_paths.size();
    
    if (VERBOSE)
      cerr << "n_nodes=" << n_nodes << endl
      << "n_sites=" << n_sites << endl;
    
    
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    std::mt19937 gen(rng_seed);
    /**************************************************************************/
    /**************************************************************************/

    // FOR NOW, WE DO NOT LEARN BRANCH LENGTHS //
      
    /**************************************************************************/
    /**************************************************************************/

    /*
    if (VERBOSE) {
      cerr << "[INITIAL GUESS AT BRANCH LENGTHS]" << endl;
      for (size_t i = 0; i < th.node_names.size(); ++i)
        cerr << th.node_names[i] << '\t' << th.branches[i] << endl;
    }
    */

    vector<double> J(8, 0.0);
    vector<double> D(8, 0.0);
    
    for (size_t i = 0; i < iteration; i++) {
      cerr << "***********************************" << endl
      << "ITR: " << i << endl;
      // estiamte J and D
      cerr << "estimating sufficient statistics" << endl;
      estimate_sufficient_statistics_by_simulation(all_paths, the_model, th,
                                                   J, D, gen);
      
      // estimate parameters
      cerr << "estimating model parameters" << endl;
      if (!OPTBRANCH) {
        compute_estimates_for_rates_only(VERBOSE, param_tol, J, D, the_model);
      }
      else
        compute_estimates_rates_and_branches(VERBOSE, param_tol, all_paths,
                                             th, the_model);
      
      estimate_root_distribution(all_paths, the_model);
      cerr << "***********************************" << endl;
    }

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    if (!out)
      throw std::runtime_error("bad output file: " + outfile);

   }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
