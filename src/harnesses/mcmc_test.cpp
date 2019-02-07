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

/* This function does the sampling for an individual change in the
   state sequence
 */
static size_t
count_jumps(const vector<vector<Path> > paths, const size_t the_site) {
  size_t counts = 0;
  for (size_t b = 1; b < paths.size(); ++b)
    counts += paths[b][the_site].jumps.size();
  return counts;
}

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
delete_branches(vector<vector<Path> > &paths, vector<string> &node_names,
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


static bool
check_homo_sites(const vector<Path> &by_site,
                 const vector<size_t> censored_site) {
  bool homo = true;
  size_t site_id = 0;
  while(homo && site_id < by_site.size()) {
    if (std::find(censored_site.begin(), censored_site.end(),
                  site_id) == censored_site.end()) {
      //homo = ((by_site[site_id].jumps.empty()) &&
      //        (!by_site[site_id].init_state)) ? true : false;
      homo = (by_site[site_id].jumps.empty()) ? true : false;
      //cerr << "checked site: " << site_id << endl;
    }
    //cerr << site_id << ": " << by_site[site_id].init_state << " "
    // << by_site[site_id].jumps.size() << " = " << homo << endl;
    site_id++;
  }
  //cerr << endl;
  return homo;
}


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


int main(int argc, const char **argv) {
  try {
    bool VERBOSE = false;
    bool FIX_NEIGHBORS = false;
    bool FIX_ROOT = false;
    bool EST_PARAM = false;
    
    string outfile;
    string statfile;
    string tracefile;

    size_t rounds = 10;
    size_t rng_seed = std::numeric_limits<size_t>::max();
    size_t the_branch = std::numeric_limits<size_t>::max();
    vector<size_t> consored_site = {2};
    size_t batch = 100;
    size_t burning = 1000;
    /* proposal: 0 - direct, 1 - unif, 2 - poisson, 3 - forward */
    size_t proposal = 0;
    
    const size_t iteration = 10; // MLE iterations
    static const double param_tol = 1e-10;
    ///////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]),
                           "collect summary stats and estimate parameters"
                           "using mcmc procedure",
                           "<param> <treefile> <path_file>");
    opt_parse.add_opt("rounds", 'r', "number of MCMC rounds",
                      false, rounds);
    opt_parse.add_opt("branch", 'b', "branch to simulate",
                      false, the_branch);
    opt_parse.add_opt("batch", 'B', "batch size",
                      false, batch);
    opt_parse.add_opt("burning", 'L', "burining length",
                      false, burning);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);
    opt_parse.add_opt("statfile", 'S', "MCMC stats",
                      false, statfile);
    opt_parse.add_opt("tracefile", 't', "MCMC trace",
                      false, tracefile);
    opt_parse.add_opt("outfile", 'o', "outfile (sampling paths)",
                      false, outfile);
    opt_parse.add_opt("estparam", 'm', "estimate model parameters",
                      false, EST_PARAM);
    opt_parse.add_opt("fix", 'f', "fix neighbors",
                      false, FIX_NEIGHBORS);
    opt_parse.add_opt("fixroot", 'R', "fix root states",
                      false, FIX_ROOT);
    opt_parse.add_opt("proposal", 'P', "MCMC proposal",
                      false, proposal);
    
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
    if (the_branch == 1) {
      // remove unwanted branches
      delete_branches(paths, node_names, th);
    }
    
    /*
    cerr << "[ROOT SEQUENCE: ]" << endl;
    for (size_t pos = 0; pos < n_sites; pos++) {
      cerr << paths[1][pos].init_state;
    }
    cerr << endl;
    */
    
    /*
    if (the_branch == 1) {
      cerr << "[LEAF SEQUENCE: ]" << endl;
      for (size_t pos = 0; pos < n_sites; pos++) {
        cerr << paths[the_branch][pos].end_state();
      }
      cerr << endl;
    }
    */
     
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
    std::ofstream of_root;
    of_root.open(rootfile.c_str());
    std::ostream out_root(of_root.rdbuf());
    
    // FILE RECORDING SUMMARY STATS
    std::ofstream outstat(statfile.c_str());
    outstat << "J_000\tJ_001\tJ_010\tJ_011\tJ_100\tJ_101\tJ_110\tJ_111\t"
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
    /* PREPARE MCMC */
    size_t itr = 0;
    size_t last_sample_point = itr;
    size_t num_sample_collected = 0;
    vector<double> pos_num_update(n_sites, 0.0);
    vector<double> J;
    vector<double> D;
    vector<vector<double> > J_batch(8, vector<double> (batch, 0.0));
    vector<vector<double> > D_batch(8, vector<double> (batch, 0.0));
    vector<double> J_accu(8, 0.0);
    vector<double> D_accu(8, 0.0);
    
    const size_t offset = 1;
    while (num_sample_collected < rounds) {

      bool update_happens = false;
      
      /* MCMC PROCEDURE */
      for (size_t site_id = offset; site_id < n_sites - offset; ++site_id) {
        vector<Path> proposed_path;
        const bool accpeted = Metropolis_Hastings_site(the_model, th, site_id,
                                                       paths, gen,
                                                       proposed_path, proposal,
                                                       FIX_ROOT);
        if (accpeted) {
          pos_num_update[site_id]++;
          update_happens = accpeted;
        }
      }
      itr++;

      /* PROCESS MCMC STATS */
      
      // CALCULATE SUFFICIENT STATS
      get_sufficient_statistics(paths, J, D);
      
      // RECORD BATCH STATS
      for (size_t i = 0; i < 8; i++) {
        J_batch[i][itr - last_sample_point - 1] = J[i];
        D_batch[i][itr - last_sample_point - 1] = D[i];
      }
      
      // NEW SAMPLING POINT
      if ((itr - last_sample_point) == batch) {
        
        if (itr > burning)
          num_sample_collected++;
        
        // COLLECT/OUTPUT J AND D
        if (itr > burning) {
          for (size_t i = 0; i < 8; i++)
            outstat << J[i] << "\t";
          for (size_t i = 0; i < 7; i++)
            outstat << D[i] << "\t";
          outstat << D[7] << endl;
          
          for(size_t i = 0; i < 8; i++) {
            J_accu[i] += J[i];
            D_accu[i] += D[i];
          }
          
          // OUTPUT ROOT
          for (size_t i = 0; i < n_sites - 1; i++)
            out_root << paths[1][i].init_state << "\t";
          out_root << paths[1][n_sites - 1].init_state << endl;
        }

        // RESET SAMPLING POINT
        last_sample_point = itr;
        
        // CALCULATE/OUTPUT BATCH STATS
        for (size_t i = 0; i < 8; i++) {
          double mean, var;
          mean_var(J_batch[i], mean, var);
          out_trace_mean << mean << "\t";
          out_trace_var << var << "\t";
        }
        for (size_t i = 0; i < 7; i++) {
          double mean, var;
          mean_var(D_batch[i], mean, var);
          out_trace_mean << mean << "\t";
          out_trace_var << var << "\t";
        }
        double mean, var;
        mean_var(D_batch[7], mean, var);
        out_trace_mean << mean << endl;
        out_trace_var << var << endl;
        
        /* (7) PARAMETER ESTIMATION */
        if (EST_PARAM) {
          EpiEvoModel updated_model = the_model;
          vector<double> J_mean(8, 0.0);
          vector<double> D_mean(8, 0.0);
          
          for (size_t i = 0; i < 8; i++) {
            J_mean[i] = J_accu[i] / batch;
            D_mean[i] = D_accu[i] / batch;
          }
          
          for (size_t i = 0; i < iteration; i++) {
            compute_estimates_for_rates_only(VERBOSE, param_tol,
                                             J_mean, D_mean, updated_model);
            estimate_root_distribution(paths, updated_model);
          }
          outparam << updated_model.T[0][0] << "\t"
          << updated_model.T[1][1] << "\t"
          << updated_model.stationary_logbaseline[0][0] << "\t"
          << updated_model.stationary_logbaseline[1][1] << "\t"
          << updated_model.init_T[0][0] << "\t"
          << updated_model.init_T[1][1] << endl;
          
          the_model = updated_model;
        }
        
        // RESET
        for (size_t i = 0; i < 8; i++) {
          J_accu[i] = 0.0;
          D_accu[i] = 0.0;
        }
      }
    }
    
    if (VERBOSE) {
      cerr << "METROPOLIS-HASTINGS NUMBER OF SUCCESS: ";
      for (size_t i = 0; i < n_sites; i++)
        cerr << pos_num_update[i] << ", ";
    }
    cerr << endl;

    /* (6) OUTPUT */
    write_root_to_pathfile_local(outfile, th.node_names.front());
    for (size_t node_id = 1; node_id < th.n_nodes; ++node_id)
      append_to_pathfile_local(outfile, th.node_names[node_id], paths[node_id]);
    
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
