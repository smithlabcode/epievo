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

/* generate initial paths by heuristics */
template <class T>
static void
initialize_paths(std::mt19937 &gen, const TreeHelper &th,
                 vector<vector<T> > &state_sequences,
                 vector<vector<Path> > &paths) {
  
  const size_t n_sites = state_sequences.front().size();
  
  paths.resize(th.n_nodes);
  
  auto unif =
  bind(std::uniform_real_distribution<double>(0.0, 1.0), std::ref(gen));
  
  vector<T> child_states = {0, 0}; // two child states
  
  for (size_t i = th.n_nodes; i > 0; --i) {
    
    const size_t node_id = i - 1;
    paths[node_id].resize(n_sites);
    
    for (size_t site_id = 0; site_id < n_sites; ++site_id) {
      
      if (!th.is_leaf(node_id)) {
        // get child states
        size_t n_ch = 0;
        for (ChildSet c(th.subtree_sizes, node_id); c.good(); ++c)
          child_states[n_ch++] = state_sequences[*c][site_id];
        
        // sample parent state
        state_sequences[node_id][site_id] =
        child_states[std::floor(unif()*n_ch)];
        
        // assign site-specific paths above two child nodes
        for (ChildSet c(th.subtree_sizes, node_id); c.good(); ++c) {
          
          const double len = th.branches[*c];
          paths[*c][site_id].tot_time = len;
          paths[*c][site_id].init_state = state_sequences[node_id][site_id];
          
          if (state_sequences[*c][site_id] != paths[*c][site_id].init_state)
            paths[*c][site_id].jumps.push_back(unif() * len);
        }
      }
    }
  }
}


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

    string outfile;               // sampled local path
    
    size_t iteration = 100;       // MCMC sample size
    size_t burnin = 10;           // burn-in MCMC iterations

    double evolutionary_time = numeric_limits<double>::lowest();

    size_t rng_seed = std::numeric_limits<size_t>::max();

    static const double param_tol = 1e-10;
    ///////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]),
                           "Simulate paths between a pair of sequences",
                           "<param> <states>");
    opt_parse.add_opt("iteration", 'i',
                      "number of MCMC-EM iteration (default: 100)",
                      false, iteration);
    opt_parse.add_opt("burnin", 'L',
                      "MCMC burn-in length (default: 10)",
                      false, burnin);
    opt_parse.add_opt("evo-time", 'T', "evolutionary time", false,
                      evolutionary_time);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);
    opt_parse.add_opt("outfile", 'o',
                      "output file of local paths (default: stdout)",
                      false, outfile);
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
    if (leftover_args.size() < 2) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    const string param_file(leftover_args[0]);
    const string states_file(leftover_args[1]);
    ///////////////////////////////////////////////////////////////////////////
    /* (1) READING PARAMETERS FROM FILE */
    if (VERBOSE)
      cerr << "[READING PARAMETERS: " << param_file << endl;
    EpiEvoModel the_model;
    read_model(param_file, the_model);
    if (VERBOSE)
      cerr << the_model << endl;

    /* (2) GET THE ROOT SEQUENCE */
    if (VERBOSE)
      cerr << "[OBTAINING ROOT/LEAF SEQUENCE]" << endl;
    TreeHelper th;
    vector<vector<bool> > state_sequences;
    read_states_file(states_file, state_sequences, th);
     
    /* (3) INITIALIZING THE RANDOM NUMBER GENERATOR */
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    std::mt19937 gen(rng_seed);
    
    /* (4) INITIALIZING THE PATH */
    vector<vector<Path> > paths;
    initialize_paths(gen, th, state_sequences, paths);
    size_t n_sites = paths.back().size();

    
    /* (5) MCMC */
    /* MCMC DATA*/
    vector<Path> proposed_path;
    // Smooth paths segmented in equally spaced intervals
    vector<double> average_path;

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    /////// STARTING MCMC
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    /* METROPOLIS-HASTINGS ALGORITHM */
    // Burning
    for(size_t burnin_itr = 0; burnin_itr < burnin; burnin_itr++) {
      for (size_t site_id = 1; site_id < n_sites - 1; ++site_id) {
        Metropolis_Hastings_site(the_model, th, site_id, paths, gen,
                                 proposed_path);
      }
    }
    
    // MCMC samples
    for(size_t mcmc_itr = 0; mcmc_itr < iteration; mcmc_itr++) {
      for (size_t site_id = 1; site_id < n_sites - 1; ++site_id) {
        Metropolis_Hastings_site(the_model, th, site_id, paths, gen,
                                 proposed_path);
        update_average_path(average_path, paths);
      }
    }
    
    /* (6) OUTPUT */
    if (VERBOSE)
      cerr << "[WRITING PATHS]" << endl;
    write_root_to_pathfile_local(outfile, th.node_names.front());
    for (size_t node_id = 1; node_id < th.n_nodes; ++node_id)
      append_to_pathfile_local(outfile, th.node_names[node_id], paths[node_id]);
  }

  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
