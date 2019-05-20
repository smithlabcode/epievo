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
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream outpath(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  if (!outpath)
    throw std::runtime_error("bad output file: " + outfile);
  
  outpath << "NODE:" << root_name << endl;
}


static void
append_to_pathfile_local(const string &pathfile, const string &node_name,
                         const vector<Path> &path_by_site) {
  std::ofstream of;
  if (!pathfile.empty()) of.open(pathfile.c_str(), std::ofstream::app);
  std::ostream outpath(pathfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  if (!outpath)
    throw std::runtime_error("bad output file: " + pathfile);
  
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
    string param_file_updated;
    string treefile_updated;

    size_t iteration = 10;           // MCMC-EM iterations
    size_t batch = 10;              // MCMC iterations
    size_t burnin = 10;           // burn-in MCMC iterations

    size_t rng_seed = std::numeric_limits<size_t>::max();

    static const double param_tol = 1e-10;
    ///////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]),
                           "Estimate parameters and evolution histories"
                           " with MCMC-EM procedure",
                           "<param> <treefile> <path_file>");
    opt_parse.add_opt("iteration", 'i',
                      "number of MCMC-EM iteration (default: 10)",
                      false, iteration);
    opt_parse.add_opt("batch", 'B',
                      "number of MCMC iteration (default: 10)",
                      false, batch);
    opt_parse.add_opt("burnin", 'L',
                      "MCMC burn-in length (default: 10)",
                      false, burnin);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);
    opt_parse.add_opt("outfile", 'o',
                      "output file of local paths (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("outparam", 'p',
                      "output file of parameters (default: stdout)",
                      false, param_file_updated);
    opt_parse.add_opt("outtree", 't',
                      "output file of tree (default: stdout)",
                      false, treefile_updated);
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
    vector<double> J;
    vector<double> D;
    vector<Path> proposed_path;
    
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    /////// START MCMC-EM
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    if (VERBOSE)
      cerr << "[RUNNING MCMC-EM (seed=" << rng_seed << ")]\n"
      << "0\t" << the_model.T[0][0] << "\t"
      << the_model.T[1][1] << "\t"
      << the_model.stationary_logbaseline[0][0] << "\t"
      << the_model.stationary_logbaseline[1][1] << "\t"
      << the_model.init_T[0][0] << "\t" << the_model.init_T[1][1] << "\t"
      << the_tree << endl;
    
    size_t current_batch = batch;
    
    /* METROPOLIS-HASTINGS ALGORITHM */
    for (size_t itr = 0; itr  < iteration; itr++) {
      // Burning
      for(size_t burnin_itr = 0; burnin_itr < burnin; burnin_itr++) {
        for (size_t site_id = 1; site_id < n_sites - 1; ++site_id) {
          Metropolis_Hastings_site(the_model, th, site_id, paths, gen,
                                   proposed_path);
        }
      }
      
      // MCMC samples
      current_batch += (itr+1) / 5;
      vector<vector<double> > J_batch(8, vector<double> (current_batch, 0.0));
      vector<vector<double> > D_batch(8, vector<double> (current_batch, 0.0));
      for(size_t mcmc_itr = 0; mcmc_itr < current_batch; mcmc_itr++) {
        for (size_t site_id = 1; site_id < n_sites - 1; ++site_id)
          Metropolis_Hastings_site(the_model, th, site_id, paths, gen,
                                   proposed_path);
        
        /* CALCULATE SUFFICIENT STATS */
        get_sufficient_statistics(paths, J, D);
        
        /* RECORD STATS (POST-BURN-IN) */
        for (size_t i = 0; i < 8; i++) {
          J_batch[i][mcmc_itr] = J[i];
          D_batch[i][mcmc_itr] = D[i];
        }
      }
      
      //////////////////////////////////////////////////////////////////////////
      ///////////               POST-MCMC PROCESSING                 ///////////
      //////////////////////////////////////////////////////////////////////////
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
      
      if (VERBOSE)
        cerr << itr+1 << "\t" << the_model.T[0][0] << "\t"
        << the_model.T[1][1] << "\t"
        << the_model.stationary_logbaseline[0][0] << "\t"
        << the_model.stationary_logbaseline[1][1] << "\t"
        << the_model.init_T[0][0] << "\t" << the_model.init_T[1][1] << "\t"
        << the_tree << endl;
    }
  
 
    /* (6) OUTPUT */
    if (VERBOSE)
      cerr << "[WRITING PATHS]" << endl;
    write_root_to_pathfile_local(outfile, th.node_names.front());
    for (size_t node_id = 1; node_id < th.n_nodes; ++node_id)
      append_to_pathfile_local(outfile, th.node_names[node_id], paths[node_id]);
    
    if (VERBOSE)
      cerr << "[WRITING PARAMETERS]\n" <<  the_model << endl;
    std::ofstream of_param;
    if (!param_file_updated.empty()) of_param.open(param_file_updated.c_str());
    std::ostream out_param(param_file_updated.empty() ?
                           std::cout.rdbuf() : of_param.rdbuf());
    if (!out_param)
      throw std::runtime_error("bad output param file: " + param_file_updated);
    out_param << the_model.format_for_param_file() << endl;
    
    if (OPTBRANCH) {
      std::ofstream of_tree;
      if (!treefile_updated.empty()) of_tree.open(treefile_updated.c_str());
      std::ostream out_tree(treefile_updated.empty() ?
                            std::cout.rdbuf() : of_tree.rdbuf());
      if (!out_tree)
        throw std::runtime_error("bad output param file: " + treefile_updated);
      out_tree << the_tree << endl;
    }

  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
