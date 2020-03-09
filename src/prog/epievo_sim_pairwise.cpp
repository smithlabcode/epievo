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
#include <sstream>
#include <algorithm>
#include <bitset>
#include <random>
#include <functional>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "Path.hpp"
#include "EpiEvoModel.hpp"
#include "EndCondSampling.hpp"
#include "ContinuousTimeMarkovModel.hpp"
#include "TreeHelper.hpp"
#include "SingleSiteSampler.hpp"
#include "TripletSampler.hpp"
#include "GlobalJump.hpp"
#include "ParamEstimation.hpp"

#include "epievo_utils.hpp"
#include "emission_utils.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::runtime_error;
using std::numeric_limits;
using std::to_string;
using std::function;
using std::placeholders::_1;
using std::bind;
using std::ofstream;

/* generate initial paths */
static void
initialize_paths_indep(std::mt19937 &gen, const vector<bool> &root_seq,
                       const vector<bool> &leaf_seq,
                       vector<vector<Path> > &paths,
                       const EpiEvoModel &the_model, const double tot_time) {

  const size_t n_sites = root_seq.size();
  paths.resize(n_sites);
  for (size_t site_id = 0; site_id < n_sites; ++site_id) {
    paths[site_id].resize(2);
  }

  std::uniform_real_distribution<double> unif(0.0, tot_time);

  paths[0][0].init_state = root_seq[0];
  paths[0][1].init_state = root_seq[0];
  paths[0][1].tot_time = tot_time;

  paths[n_sites-1][0].init_state = root_seq[n_sites-1];
  paths[n_sites-1][1].init_state = root_seq[n_sites-1];
  paths[n_sites-1][1].tot_time = tot_time;

  if (paths[0][1].end_state() != leaf_seq[0])
    paths[0][1].jumps.push_back(unif(gen));
  if (paths[n_sites-1][1].end_state() != leaf_seq[n_sites-1])
    paths[n_sites-1][1].jumps.push_back(unif(gen));

  for (size_t site_id = 1; site_id < n_sites - 1; ++site_id) {
    const size_t l = site_id - 1;
    const size_t r = site_id + 1;
    const size_t pattern0 = triple2idx(root_seq[l], false, root_seq[r]);
    const size_t pattern1 = triple2idx(root_seq[l], true, root_seq[r]);
    const bool start_state = root_seq[site_id];
    const bool end_state = leaf_seq[site_id];
    const TwoStateCTMarkovModel ctmm(the_model.triplet_rates[pattern0],
                                     the_model.triplet_rates[pattern1]);

    paths[site_id][0].init_state = start_state;
    paths[site_id][1].init_state = start_state;
    paths[site_id][1].tot_time = tot_time;

    end_cond_sample_forward_rejection(ctmm, start_state, end_state, tot_time,
                                      gen, paths[site_id][1].jumps, 0.0);
    assert(paths[site_id][1].end_state() == end_state);
  }
}

static void
write_root_to_pathfile_local(const string &outfile, const string &root_name) {
  ofstream out(outfile);
  if (!out)
    throw std::runtime_error("bad output file: " + outfile);
  out << "NODE:" << root_name << endl;
}

static void
append_to_pathfile_local(const string &pathfile, const string &node_name,
                         const vector<vector<Path> > &paths,
                         const size_t node_id) {
  std::ofstream out(pathfile, std::ofstream::app);
  if (!out)
    throw std::runtime_error("bad output file: " + pathfile);

  out << "NODE:" << node_name << endl;
  for (size_t i = 0; i < paths.size(); ++i)
    out << i << '\t' << paths[i][node_id] << '\n';
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char **argv) {
  try {

    bool VERBOSE = false;
    string outfile;               // sampled local path
    string pathfile;               // initial local path
    size_t burnin = 10;           // burn-in MCMC iterations
    double evolutionary_time = 1.0;

    size_t rng_seed = std::numeric_limits<size_t>::max();

    ///////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]),
                           "simulate paths between a pair of sequences",
                           "<param> <states>");
    opt_parse.add_opt("burnin", 'L', "MCMC burn-in length (default: 10)",
                      false, burnin);
    opt_parse.add_opt("evo-time", 'T', "evolutionary time", false,
                      evolutionary_time);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);
    opt_parse.add_opt("outfile", 'o', "local paths output file", true, outfile);
    opt_parse.add_opt("paths", 'p', "local paths input file", false, pathfile);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);

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
    the_model.scale_triplet_rates();
    the_model.use_init_T = false;

    if (VERBOSE)
      cerr << the_model << endl;

    /* (2) GET THE ROOT SEQUENCE */
    if (VERBOSE)
      cerr << "[OBTAINING ROOT/LEAF SEQUENCE]" << endl;
    TreeHelper th(evolutionary_time);
    vector<string> node_names; // might get used twice
    vector<state_seq> node_seqs;
    read_states_file(states_file, node_names, node_seqs);
    if (node_seqs.size() != 2)
      throw runtime_error("require exactly 2 nodes in: " + states_file);

    const state_seq root_seq(std::move(node_seqs.front()));
    const state_seq leaf_seq(std::move(node_seqs.back()));

    const size_t n_sites = root_seq.size();
    /* (3) INITIALIZING THE RANDOM NUMBER GENERATOR */
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    std::mt19937 gen(rng_seed);

    vector<vector<Path> > paths; // [sites] x [nodes]
    if (!pathfile.empty()) {
      // reading into node_names again here
      vector<vector<Path> > orig_paths; // [nodes] x [sites]
      read_paths(pathfile, node_names, orig_paths);
      
      ////////// transposing
      paths.resize(orig_paths[1].size());
      for (size_t i = 0; i < orig_paths[1].size(); ++i) {
        paths[i].resize(orig_paths.size());
        for (size_t b = 1; b < 2; ++b) {
          paths[i][b] = orig_paths[b][i];
        }
      }
      //////////
    }
    else
      initialize_paths_indep(gen, root_seq, leaf_seq, paths, the_model,
                             evolutionary_time);

    if (VERBOSE)
      cerr << "[SIMULATING]" << endl;

    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    /////// STARTING MCMC
    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    vector<vector<vector<double> > > emit(n_sites);
    for (size_t site_id = 0; site_id < n_sites; site_id++) {
      emit[site_id].resize(2);
      emit[site_id][0].resize(2);
      emit[site_id][0][0] = (paths[site_id][1].init_state ? 0.0 : 1.0);
      emit[site_id][0][1] = (paths[site_id][1].init_state ? 1.0 : 0.0);
    }

    vector<double> tri_llh(n_sites, 0.0); // log-likelihood over three triplets
    // pre-compute triplet log-likelihood on each site
    for (size_t site_id = 1; site_id < n_sites - 1; ++site_id)
      tri_llh[site_id] = path_log_likelihood(the_model, paths[site_id-1],
                                             paths[site_id], paths[site_id+1]);
    /* METROPOLIS-HASTINGS ALGORITHM */
    // Burning
    for (size_t burnin_itr = 0; burnin_itr < burnin; burnin_itr++) {
      for (size_t site_id = 1; site_id < n_sites - 1; ++site_id) {
        Metropolis_Hastings_site(the_model, th, site_id, paths,
                                 emit[site_id], tri_llh[site_id-1],
                                 tri_llh[site_id], tri_llh[site_id+1], gen);
      }
    }
    write_root_to_pathfile_local(outfile, th.node_names.front());
    append_to_pathfile_local(outfile, th.node_names[1], paths, 1);

    /* (6) OUTPUT */
    if (VERBOSE)
      cerr << "[WRITING OUTPUT]" << endl;
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
