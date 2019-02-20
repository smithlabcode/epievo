/* Copyright (C) 2019 University of Southern California
 *                    Xiaojing Ji, Jianghan Qu and Andrew D Smith
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
#include <algorithm>
#include <bitset>
#include <random>
#include <functional>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "Path.hpp"
#include "EpiEvoModel.hpp"
#include "StateSeq.hpp"
#include "EndCondSampling.hpp"
#include "ContinuousTimeMarkovModel.hpp"
#include "TreeHelper.hpp"
#include "SingleSiteSampler.hpp"
#include "TripletSampler.hpp"
#include "GlobalJump.hpp"
#include "ParamEstimation.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::runtime_error;
using std::numeric_limits;

///////////////////////////////////////////////////////////////////////////

static void
write_root_to_pathfile_local(const string &outfile, const string &root_name) {
  std::ofstream outpath(outfile.c_str());
  outpath << "NODE:" << root_name << endl;
}


static void
append_to_pathfile_local(const string &pathfile, const string &node_name,
                         const vector<Path> &path_by_site) {
  std::ofstream outpath(pathfile.c_str(), std::ofstream::app);
  outpath << "NODE:" << node_name << endl;
  for (size_t i = 0; i < path_by_site.size(); ++i)
    outpath << i << '\t' << path_by_site[i] << '\n';
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char **argv) {
  try {
    bool VERBOSE = false;
    bool OPTBRANCH = false;

    string outfile;

    size_t iteration = 10;           // MCMC-EM iterations
    size_t batch = 10;              // MCMC iterations
    size_t burnin = 10;           // burn-in MCMC iterations

    size_t rng_seed = std::numeric_limits<size_t>::max();

    static const double param_tol = 1e-10;
    ///////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]),
                           "Collect summary stats and estimate parameters"
                           "with MCMC-EM procedure",
                           "<param> <treefile> <path_file>");
    opt_parse.add_opt("iteration", 'i', "number of MCMC iteration",
                      false, iteration);
    opt_parse.add_opt("batch", 'B', "batch size", false, batch);
    opt_parse.add_opt("burnin", 'L', "burining length", false, burnin);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);
    opt_parse.add_opt("outfile", 'o', "outfile (sampling paths)",
                      false, outfile);
    opt_parse.add_opt("branch", 'b', "optimize branch lengths as well",
                      false, OPTBRANCH);
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
    if (leftover_args.size() < 3) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    const string param_file(leftover_args[0]);
    const string treefile(leftover_args[1]);
    const string input_file(leftover_args[2]);
    ///////////////////////////////////////////////////////////////////////////
    /* (1) LOADING (FAKE) TREE */
    if (VERBOSE)
      cerr << "[READING TREE: " << treefile << "]" << endl;
    PhyloTreePreorder the_tree; // tree topology and branch lengths
    std::ifstream tree_in(treefile.c_str());
    if (!tree_in || !(tree_in >> the_tree))
      throw runtime_error("cannot read tree file: " + treefile);
    TreeHelper th(the_tree);

    /* (2) READING PARAMETERS FROM FILE */
    if (VERBOSE)
      cerr << "[READING PARAMETERS: " << param_file << endl;
    EpiEvoModel the_model;
    read_model(param_file, the_model);
    if (VERBOSE)
      cerr << the_model << endl;

    /* (3) LOADING PATHS */
    if (VERBOSE)
      cerr << "[READING PATHS FILE: " << input_file << "]" << endl;
    vector<string> node_names;
    vector<vector<Path> > paths; // along multiple branches
    read_paths(input_file, node_names, paths);
    const size_t n_sites = paths[1].size();
     
    /* (4) INITIALIZING THE RANDOM NUMBER GENERATOR */
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    std::mt19937 gen(rng_seed);
    if (VERBOSE)
      cerr << "[INITIALIZED RNG (seed=" << rng_seed << ")]\n"
      << "itr\tstationary\tbaseline\tinit\ttree" << endl;
    
    /* (5) MCMC */
  
    /* MCMC PARAMETERS AND DATA*/
    size_t mcmc_itr = 0;                   // MCMC iterations elapsed
    size_t finished_iterations = 0;        // finished MCMC-EM iterations
    vector<double> J;
    vector<double> D;
    vector<vector<double> > J_batch(8, vector<double> (batch, 0.0));
    vector<vector<double> > D_batch(8, vector<double> (batch, 0.0));
    vector<Path> proposed_path;
    
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    /////// START MCMC-EM
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    if (VERBOSE)
      cerr << "[RUNNING MCMC-EM (seed=" << rng_seed << ")]\n"
      << finished_iterations << "\t" << the_model.T[0][0] << "\t"
      << the_model.T[1][1] << "\t"
      << the_model.stationary_logbaseline[0][0] << "\t"
      << the_model.stationary_logbaseline[1][1] << "\t"
      << the_model.init_T[0][0] << "\t" << the_model.init_T[1][1] << "\t"
      << the_tree << endl;
    
    while (finished_iterations < iteration) {
      
      /* METROPOLIS-HASTINGS ALGORITHM */
      for (size_t site_id = 1; site_id < n_sites - 1; ++site_id) {
        Metropolis_Hastings_site(the_model, th, site_id, paths, gen,
                                 proposed_path);
      }
      mcmc_itr++;
      //////////////////////////////////////////////////////////////////////////
      ///////////               POST-MCMC PROCESSING                 ///////////
      //////////////////////////////////////////////////////////////////////////
      /* CALCULATE SUFFICIENT STATS */
      get_sufficient_statistics(paths, J, D);
      
      if (mcmc_itr > burnin)
      /* RECORD STATS (POST-BURN-IN) */
        for (size_t i = 0; i < 8; i++) {
          J_batch[i][mcmc_itr - burnin - 1] = J[i];
          D_batch[i][mcmc_itr - burnin - 1] = D[i];
        }
      
      /* MCMC OVER */
      if (mcmc_itr == (burnin + batch)) {

        /* CALCULATE/OUTPUT BATCH AVERAGE/VAR */
        vector<double> J_mean(8, 0.0);
        vector<double> D_mean(8, 0.0);

        for (size_t i = 0; i < 8; i++) {
          J_mean[i] = std::accumulate(J_batch[i].begin(), J_batch[i].end(), 0.0)
          / J_batch[i].size();
          D_mean[i] = std::accumulate(D_batch[i].begin(), D_batch[i].end(), 0.0)
          / D_batch[i].size();
        }
        
        /* PARAMETER ESTIMATION */
        if (!OPTBRANCH)
          compute_estimates_for_rates_only(false, param_tol,
                                           J_mean, D_mean, the_model);
        else {
          compute_estimates_rates_and_branches(false, param_tol, paths,
                                               th, the_model);
          the_tree.set_branch_lengths(th.branches);
        }

        estimate_root_distribution(paths, the_model);
        
        /* RESET MCMC BATCH */
        mcmc_itr = 0;
        finished_iterations++;
        
        if (VERBOSE)
          cerr << finished_iterations << "\t" << the_model.T[0][0] << "\t"
          << the_model.T[1][1] << "\t"
          << the_model.stationary_logbaseline[0][0] << "\t"
          << the_model.stationary_logbaseline[1][1] << "\t"
          << the_model.init_T[0][0] << "\t" << the_model.init_T[1][1] << "\t"
          << the_tree << endl;
      }
    }

    /* (6) OUTPUT GLOBAL PATH */
    write_root_to_pathfile_local(outfile, th.node_names.front());
    for (size_t node_id = 1; node_id < th.n_nodes; ++node_id)
      append_to_pathfile_local(outfile, th.node_names[node_id], paths[node_id]);
    
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
