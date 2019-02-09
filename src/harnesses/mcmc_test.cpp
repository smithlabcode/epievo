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


static void
print_summary_stats(std::ofstream &out, const vector<double> &J,
                    const vector<double> &D) {
  for (size_t i = 0; i < 8; i++)
    out << J[i] << "\t";
  for (size_t i = 0; i < 7; i++)
    out << D[i] << "\t";
  out << D[7] << endl;
}

static void
print_root_states(std::ofstream &out, const vector<vector<Path> > &paths) {
  for (size_t i = 0; i < paths[1].size(); i++)
    out << paths[1][i].init_state << "\t";
  out << paths[1][paths[1].size() - 1].init_state << endl;
}


/* This function will erase all tree branches but only keep the first one.
 */
static void
delete_nonfirst_branches(vector<vector<Path> > &paths, vector<string> &node_names,
                   TreeHelper &th) {
  while(paths.size() > 2) {
    paths.pop_back();
    node_names.pop_back();
    th.branches.pop_back();
    th.parent_ids.pop_back();
    th.subtree_sizes.pop_back();
    th.subtree_sizes[0]--;
    th.n_nodes--;
  }
}


/* Compute mean and variance given a vector.
 */
static void
mean_var(const vector<double> x, double &mean, double &var) {
  mean = 0;
  var = 0;
  const double sq_sum = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);

  if (x.size() > 0) {
    mean = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
    var = sq_sum / x.size() - mean * mean;
  }
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char **argv) {
  try {
    bool VERBOSE = false;
    bool FIX_ROOT = false;
    bool EST_PARAM = false;
    bool TEST_FIRST_BRANCH = false;  // only test the first branch

    string outfile;
    string statfile;
    string tracefile;

    size_t iteration = 10;           // MCMC-EM iterations
    const size_t iteration_mle = 10; // MLE iterations
    size_t batch = 100;              // MCMC iterations
    size_t burnin = 1000;           // burn-in MCMC iterations

    size_t rng_seed = std::numeric_limits<size_t>::max();

    /* proposal: 0 - poisson, 1 - direct, 2 - forward, 3 - unif */
    size_t proposal = 0;
    
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
    opt_parse.add_opt("statfile", 'S', "summary stats", false, statfile);
    opt_parse.add_opt("tracefile", 't', "MCMC trace", false, tracefile);
    opt_parse.add_opt("outfile", 'o', "outfile (sampling paths)",
                      false, outfile);
    opt_parse.add_opt("branch", 'b', "only test the first branch",
                      false, TEST_FIRST_BRANCH);
    opt_parse.add_opt("estparam", 'm', "estimate model parameters",
                      false, EST_PARAM);
    opt_parse.add_opt("fixroot", 'R', "fix root states",
                      false, FIX_ROOT);
    opt_parse.add_opt("proposal", 'P', "MCMC proposal choice",
                      false, proposal);
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
    if (leftover_args.size() < 1) {
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
    the_model.scale_triplet_rates();
    if (VERBOSE)
      cerr << the_model << endl;

    /* (3) LOADING PATHS */
    if (VERBOSE)
      cerr << "[READING PATHS FILE: " << input_file << "]" << endl;
    vector<string> node_names;
    vector<vector<Path> > paths; // along multiple branches
    read_paths(input_file, node_names, paths);
    const size_t n_sites = paths[1].size();
    if (TEST_FIRST_BRANCH)       // if only want to test single branch
      delete_nonfirst_branches(paths, node_names, th);

     
    /* (4) INITIALIZING THE RANDOM NUMBER GENERATOR */
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    std::mt19937 gen(rng_seed);
    if (VERBOSE)
      cerr << "[INITIALIZED RNG (seed=" << rng_seed << ")]" << endl;
    
    /* (5) MCMC */
    // FILE RECORDING ROOT STATES
    string rootfile = statfile + ".root";
    std::ofstream out_root(rootfile.c_str());
    
    // FILE RECORDING SUMMARY STATS
    std::ofstream out_stat(statfile.c_str());
    out_stat << "J_000\tJ_001\tJ_010\tJ_011\tJ_100\tJ_101\tJ_110\tJ_111\t"
    << "D_000\tD_001\tD_010\tD_011\tD_100\tD_101\tD_110\tD_111" << endl;

    // FILE RECORDING MODEL PARAMETERS
    string outparamfile = statfile + ".param";;
    std::ofstream outparam(outparamfile.c_str());
    outparam << "st0\tst1\tbl0\tbl1\trt0\trt1" << endl;
    
    // FILE RECORDING BATCH MEAN/VAR
    string trace_mean_file = tracefile + ".mean";
    string trace_var_file = tracefile + ".var";
    std::ofstream of_trace_mean, of_trace_var;
    of_trace_mean.open(trace_mean_file.c_str());
    of_trace_var.open(trace_var_file.c_str());

    std::ostream out_trace_mean(of_trace_mean.rdbuf());
    std::ostream out_trace_var(of_trace_var.rdbuf());
    out_trace_mean << "J_000\tJ_001\tJ_010\tJ_011\tJ_100\tJ_101\tJ_110\tJ_111\t"
    << "D_000\tD_001\tD_010\tD_011\tD_100\tD_101\tD_110\tD_111" << endl;
    out_trace_var << "J_000\tJ_001\tJ_010\tJ_011\tJ_100\tJ_101\tJ_110\tJ_111\t"
    << "D_000\tD_001\tD_010\tD_011\tD_100\tD_101\tD_110\tD_111" << endl;

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

    while (finished_iterations < iteration) {
      
      /* METROPOLIS-HASTINGS ALGORITHM */
      for (size_t site_id = 1; site_id < n_sites - 1; ++site_id) {
        // XJ: We current turn off the interface of choosing Metropolis-Hastings
        //     proposal. Poisson proposal is always used.
        Metropolis_Hastings_site(the_model, th, site_id, paths, gen,
                                 proposed_path, FIX_ROOT);
      }
      mcmc_itr++;

      //////////////////////////////////////////////////////////////////////////
      /////// PROCESS MCMC STATISTICS
      //////////////////////////////////////////////////////////////////////////

      /* CALCULATE SUFFICIENT STATS */
      get_sufficient_statistics(paths, J, D);
      
      if (mcmc_itr > burnin)
      /* RECORD STATS (POST-BURN-IN) */
        for (size_t i = 0; i < 8; i++) {
          J_batch[i][mcmc_itr - burnin - 1] = J[i];
          D_batch[i][mcmc_itr - burnin - 1] = D[i];
        }
      
      /* FINISH MCMC BATCH */
      if (mcmc_itr == (burnin + batch)) {
        /* PRINT STATS */
        print_summary_stats(out_stat, J, D);
        print_root_states(out_root, paths);

        /* CALCULATE/OUTPUT BATCH AVERAGE/VAR */
        vector<double> J_mean(8, 0.0);
        vector<double> D_mean(8, 0.0);
        vector<double> J_var(8, 0.0);
        vector<double> D_var(8, 0.0);
        for (size_t i = 0; i < 8; i++) {
          mean_var(J_batch[i], J_mean[i], J_var[i]);
          out_trace_mean << J_mean[i] << "\t";
          out_trace_var << J_var[i] << "\t";
        }
        for (size_t i = 0; i < 7; i++) {
          mean_var(D_batch[i], D_mean[i], D_var[i]);
          out_trace_mean << D_mean[i] << "\t";
          out_trace_var << D_var[i] << "\t";
        }
        mean_var(D_batch[7], D_mean[7], D_var[7]);
        out_trace_mean << D_mean[7] << endl;
        out_trace_var << D_var[7] << endl;
        
        /* PARAMETER ESTIMATION */
        if (EST_PARAM) {
          for (size_t i = 0; i < iteration_mle; i++) {
            compute_estimates_for_rates_only(VERBOSE, param_tol,
                                             J_mean, D_mean, the_model);
            estimate_root_distribution(paths, the_model);
          }
          outparam << the_model.T[0][0] << "\t"
          << the_model.T[1][1] << "\t"
          << the_model.stationary_logbaseline[0][0] << "\t"
          << the_model.stationary_logbaseline[1][1] << "\t"
          << the_model.init_T[0][0] << "\t"
          << the_model.init_T[1][1] << endl;
        }
        
        /* RESET BATCH */
        mcmc_itr = 0;
        finished_iterations++;
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
