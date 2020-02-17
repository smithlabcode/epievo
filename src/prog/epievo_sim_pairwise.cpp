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
#include "StateSeq.hpp"
#include "EpiEvoModel.hpp"
#include "StateSeq.hpp"
#include "EndCondSampling.hpp"
#include "ContinuousTimeMarkovModel.hpp"
#include "TreeHelper.hpp"
#include "SingleSiteSampler.hpp"
#include "EmitDistro.hpp"
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
using std::to_string;
using std::function;
using std::placeholders::_1;
using std::bind;

///////////////////////////////////////////////////////////////////////////
template <class T>
size_t
read_states_file(const string &statesfile, vector<T> &root_seq,
                 vector<T> &leaf_seq) {

  std::ifstream in(statesfile.c_str());
  if (!in)
    throw std::runtime_error("bad states file: " + statesfile);

  // first line is the list of node names
  string buffer;
  if (!getline(in, buffer))
    throw std::runtime_error("cannot read nodes line in: " + statesfile);

  std::istringstream nodes_iss(buffer);
  vector<string> node_names;

  string tmp_node_name;
  while ((node_names.size() < 2) && (nodes_iss >> tmp_node_name))
    node_names.push_back(tmp_node_name);

  if (node_names.size() != 2)
    throw std::runtime_error("fewer than 2 nodes in: " + statesfile);

  if (node_names.front()[0] == '#')
    node_names.front() = node_names.front().substr(1);

  // read states
  size_t site_count = 0;
  while (getline(in, buffer)) {
    std::istringstream iss;
    iss.str(std::move(buffer));

    size_t site_index = 0;
    iss >> site_index; // not important info but must be removed

    // now read the states for the current site
    T tmp_state_val;
    iss >> tmp_state_val; // the first state is root
    root_seq.push_back(tmp_state_val);
    iss >> tmp_state_val; // the second state is leaf
    leaf_seq.push_back(tmp_state_val);
    ++site_count;
  }

  if (site_count == 0)
    throw std::runtime_error("no sites read from states file: " + statesfile);

  return site_count;
}


/* generate initial paths */
static void
initialize_paths_indep(std::mt19937 &gen, const vector<bool> &root_seq,
                       const vector<bool> &leaf_seq,
                       vector<vector<Path> > &paths,
                       const EpiEvoModel &the_model, const double tot_time) {

  const size_t n_sites = root_seq.size();
  paths.resize(2);
  paths[0].resize(n_sites);
  paths[1].resize(n_sites);

  std::uniform_real_distribution<double> unif(0.0, tot_time);

  paths[0][0].init_state = root_seq[0];
  paths[1][0].init_state = root_seq[0];
  paths[1][0].tot_time = tot_time;

  paths[0][n_sites-1].init_state = root_seq[n_sites-1];
  paths[1][n_sites-1].init_state = root_seq[n_sites-1];
  paths[1][n_sites-1].tot_time = tot_time;

  if (paths[1][0].end_state() != leaf_seq[0])
    paths[1][0].jumps.push_back(unif(gen));
  if (paths[1][n_sites-1].end_state() != leaf_seq[n_sites-1])
    paths[1][n_sites-1].jumps.push_back(unif(gen));

  for (size_t site_id = 1; site_id < n_sites - 1; ++site_id) {
    const size_t pattern0 = triple2idx(root_seq[site_id-1], false,
                                       root_seq[site_id+1]);
    const size_t pattern1 = triple2idx(root_seq[site_id-1], true,
                                       root_seq[site_id+1]);
    const bool start_state = root_seq[site_id];
    const bool end_state = leaf_seq[site_id];
    const TwoStateCTMarkovModel ctmm(the_model.triplet_rates[pattern0],
                                     the_model.triplet_rates[pattern1]);

    paths[0][site_id].init_state = start_state;
    paths[1][site_id].init_state = start_state;
    paths[1][site_id].tot_time = tot_time;

    end_cond_sample_forward_rejection(ctmm, start_state, end_state, tot_time,
                                      gen, paths[1][site_id].jumps, 0.0);
    assert(paths[1][site_id].end_state() == end_state);
  }
}


static void
write_root_to_pathfile_local(const string &outfile, const string &root_name) {
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile);
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
                           "Simulate paths between a pair of sequences",
                           "<param> <states>");
    opt_parse.add_opt("burnin", 'L',
                      "MCMC burn-in length (default: 10)",
                      false, burnin);
    opt_parse.add_opt("evo-time", 'T', "evolutionary time", false,
                      evolutionary_time);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);
    opt_parse.add_opt("outfile", 'o',
                      "output file of local paths (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("pathfile", 'p',
                      "input file of local paths", false, pathfile);
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
    the_model.scale_triplet_rates();
    the_model.use_init_T = false;

    if (VERBOSE)
      cerr << the_model << endl;

    /* (2) GET THE ROOT SEQUENCE */
    if (VERBOSE)
      cerr << "[OBTAINING ROOT/LEAF SEQUENCE]" << endl;
    TreeHelper th(evolutionary_time);
    vector<bool> root_seq, leaf_seq;
    read_states_file(states_file, root_seq, leaf_seq);
    const size_t n_sites = root_seq.size();
    /* (3) INITIALIZING THE RANDOM NUMBER GENERATOR */
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    std::mt19937 gen(rng_seed);

    vector<vector<Path> > paths;
    vector<string> node_names;
    if (!pathfile.empty())
      read_paths(pathfile, node_names, paths);
    else
      initialize_paths_indep(gen, root_seq, leaf_seq, paths, the_model,
                             evolutionary_time);
    if (paths.size() > 2) // truncate extra nodes
      paths.resize(2);

    if (VERBOSE)
      cerr << "[SIMULATING]" << endl;

    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    /////// STARTING MCMC
    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    /*
    vector<vector<vector<double> > > all_emit; // [sites x nodes x 2]
    const vector<function<double (bool)> > emit_root = {Bernoulli(0.1),
      Bernoulli(0.95)};
    const vector<function<double (bool)> > emit_leaf = {Bernoulli(0.1),
      Bernoulli(0.95)};
    compute_emit(root_seq, leaf_seq, all_emit, emit_root, emit_leaf);
    */
    vector<vector<vector<double> > > emit(n_sites);
    for (size_t site_id = 0; site_id < n_sites; site_id++) {
      emit[site_id].resize(2);
      emit[site_id][0].resize(2);
      emit[site_id][0][0] = (paths[1][site_id].init_state == false ?
                             1.0 : 0.0);
      emit[site_id][0][1] = (paths[1][site_id].init_state == true ?
                             1.0 : 0.0);
    }

    /* METROPOLIS-HASTINGS ALGORITHM */
    // Burning
    for(size_t burnin_itr = 0; burnin_itr < burnin; burnin_itr++) {
      for (size_t site_id = 1; site_id < n_sites - 1; ++site_id) {
        Metropolis_Hastings_site(the_model, th, site_id, paths,
                                 emit[site_id], gen);
      }
    }
    write_root_to_pathfile_local(outfile, th.node_names.front());
    append_to_pathfile_local(outfile, th.node_names[1], paths[1]);

    /* (6) OUTPUT */
    if (VERBOSE)
      cerr << "[WRITING OUTPUT]" << endl;
  }

  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
