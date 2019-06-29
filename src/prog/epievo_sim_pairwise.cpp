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
#include "IntervalSampler.hpp"
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

///////////////////////////////////////////////////////////////////////////
template <class T>
size_t
read_states_file(const string &statesfile, vector<vector<T> > &state_sequences,
                 TreeHelper &th) {
  
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
  while (nodes_iss >> tmp_node_name)
    node_names.push_back(tmp_node_name);
  
  if (node_names.size() != 2)
    throw std::runtime_error("more or fewer than 2 nodes in: " + statesfile);
  
  if (node_names.front()[0] == '#') {
    if (node_names.front().length() == 1)
      node_names = vector<string>(node_names.begin() + 1, node_names.end());
    else node_names.front() = node_names.front().substr(1);
  }
  
  th.node_names = node_names;
  
  // read states
  const size_t n_nodes = node_names.size();
  size_t site_count = 0;
  
  state_sequences.resize(th.n_nodes);
  
  while (getline(in, buffer)) {
    std::istringstream iss;
    iss.str(std::move(buffer));
    
    size_t site_index = 0;
    iss >> site_index; // not important info but must be removed
    
    // now read the states for the current site
    size_t node_idx = 0;
    T tmp_state_val;
    while (node_idx < n_nodes && iss >> tmp_state_val)
        state_sequences[node_idx++].push_back(tmp_state_val);
    
    if (node_idx < n_nodes)
      throw std::runtime_error("inconsistent number of states: " +
                               to_string(node_idx) + "/" + to_string(n_nodes));
    ++site_count;
  }
  
  if (site_count == 0)
    throw std::runtime_error("no sites read from states file: " + statesfile);

  return site_count;
}



/* generate initial paths */
static void
initialize_paths_indep(std::mt19937 &gen,
                       const vector<vector<size_t> > &state_sequences,
                       vector<vector<Path> > &paths,
                       const EpiEvoModel &the_model, const double tot_time) {

  const size_t n_sites = state_sequences.front().size();
  paths.resize(2);
  paths[0].resize(n_sites);
  paths[1].resize(n_sites);
  
  const vector<double> triplet_rates = the_model.triplet_rates;

  std::uniform_real_distribution<double> unif(0.0, tot_time);
  paths[0][0].init_state = state_sequences[0][0];
  paths[1][0].init_state = state_sequences[0][0];
  paths[1][0].tot_time = tot_time;
  paths[0][n_sites-1].init_state = state_sequences[0][n_sites-1];
  paths[1][n_sites-1].init_state = state_sequences[0][n_sites-1];
  paths[1][n_sites-1].tot_time = tot_time;
  if (paths[1][0].end_state() != state_sequences[1][0])
    paths[1][0].jumps.push_back(unif(gen));
  if (paths[1][n_sites-1].end_state() != state_sequences[1][n_sites-1])
    paths[1][n_sites-1].jumps.push_back(unif(gen));
  
  for (size_t site_id = 1; site_id < n_sites - 1; ++site_id) {
    const size_t pattern0 = triple2idx(state_sequences[0][site_id-1], false,
                                        state_sequences[0][site_id+1]);
    const size_t pattern1 = triple2idx(state_sequences[0][site_id-1], true,
                                        state_sequences[0][site_id+1]);
    const size_t start_state = state_sequences[0][site_id];
    const size_t end_state = state_sequences[1][site_id];
    // cerr << triplet_rates[pattern0] << ", " << triplet_rates[pattern1] << endl; 
    const TwoStateCTMarkovModel ctmm(triplet_rates[pattern0], triplet_rates[pattern1]);
    
    paths[0][site_id].init_state = start_state;
    paths[1][site_id].init_state = start_state;
    paths[1][site_id].tot_time = tot_time;
    
    //if (paths[1][site_id].end_state() != end_state) {
      //std::uniform_real_distribution<double> unif(0.0, tot_time);
      //const double u = unif(gen);
    //  paths[1][site_id].jumps.push_back(tot_time/2);
    //}
    end_cond_sample_Poisson(ctmm, start_state, end_state,
                            tot_time, gen, paths[1][site_id].jumps, 0.0);
    assert(paths[1][site_id].end_state() == end_state);
  }
}


static void
add_paths(const vector<vector<Path> > &paths,
          vector<vector<double> > &average_paths) {
  const size_t n_points = average_paths.front().size();
  const double bin = paths[1][0].tot_time / (n_points - 1);
  for (size_t site_id = 0; site_id < paths.back().size(); site_id++) {
    size_t prev_state = paths[1][site_id].init_state;
    average_paths[site_id][0] += prev_state;
    double curr_time = bin;
    for (size_t i = 1; i < average_paths[site_id].size(); i++) {
        average_paths[site_id][i] +=
      paths[1][site_id].state_at_time(curr_time);
        curr_time += bin;
    }
  }
}


template<class T>
double vector_distance(vector<T> v1, vector<T> v2) {
  double ret = 0.0;
  for(size_t i = 0; i < v1.size(); i++) {
    const double dist = (v1[i] - v2[i]);
    ret += dist * dist;
  }
  return ret > 0.0 ? sqrt(ret) : 0.0;
}

      
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

int main(int argc, const char **argv) {
  try {
    bool VERBOSE = false;

    string outfile;               // sampled local path
    string pathfile;               // initial local path
    
    size_t n_samples = 100;         // number of paths to simulate
    size_t burnin = 10;           // burn-in MCMC iterations
    // number of time points to output (including two ending points)
    size_t n_points = 100;
    
    double evolutionary_time = 1.0;

    size_t rng_seed = std::numeric_limits<size_t>::max();

    ///////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]),
                           "Simulate paths between a pair of sequences",
                           "<param> <states>");
    opt_parse.add_opt("number", 'i',
                      "number of paths to simulate (default: 100)",
                      false, n_samples);
    opt_parse.add_opt("points", 'n',
                      "number of time points to output (default: 100)",
                      false, n_points);
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

    if (VERBOSE)
      cerr << the_model << endl;

    /* (2) GET THE ROOT SEQUENCE */
    if (VERBOSE)
      cerr << "[OBTAINING ROOT/LEAF SEQUENCE]" << endl;
    TreeHelper th(evolutionary_time);
    vector<vector<size_t> > state_sequences; // double if data is continuous
    read_states_file(states_file, state_sequences, th);
    const size_t n_sites = state_sequences.front().size();
    
    /* (3) INITIALIZING THE RANDOM NUMBER GENERATOR */
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    std::mt19937 gen(rng_seed);
    
    vector<vector<double> > average_paths =
    vector<vector<double> > (n_sites, vector<double> (n_points, 0.0));

    vector<vector<Path> > init_paths;
    vector<string> node_names;
    if (!pathfile.empty())
      read_paths(pathfile, node_names, init_paths);
    else 
      initialize_paths_indep(gen, state_sequences, init_paths, the_model,
                             evolutionary_time);

    if (VERBOSE)
      cerr << "[SIMULATING]" << endl;
    if (n_samples == 0) {
      vector<vector<Path> > paths = init_paths;
      add_paths(paths, average_paths);
    }

    for(size_t sample_id = 0; sample_id < n_samples; sample_id++) {
      if (VERBOSE)
        cout << "sample: " << sample_id;
      
      /* INITIALIZING THE PATH */
      vector<vector<Path> > paths = init_paths;
      size_t n_sites = paths.back().size();
      //////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
      /////// STARTING MCMC
      //////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////

      /* METROPOLIS-HASTINGS ALGORITHM */
      // Burning
      double n_acc = 0;
      for(size_t burnin_itr = 0; burnin_itr < burnin; burnin_itr++) {
        for (size_t site_id = 1; site_id < n_sites - 1; ++site_id) {
          vector<Path> proposed_path;
          //for(size_t i = 1; i < paths[1][site_id-1].jumps.size(); ++i)
          //  assert(paths[1][site_id-1].jumps[i] > paths[1][site_id-1].jumps[i-1]);
          //for(size_t i = 1; i < paths[1][site_id+1].jumps.size(); ++i)
          //  assert(paths[1][site_id+1].jumps[i] > paths[1][site_id+1].jumps[i-1]);
          n_acc += Metropolis_Hastings_interval(the_model, th, site_id, paths,
                                                gen);
          assert(paths[1][site_id].init_state == state_sequences[0][site_id]);
          assert(paths[1][site_id].end_state() == state_sequences[1][site_id]);
        }
        //vector<double> J_new, D_new;
        //get_sufficient_statistics(paths, J_new, D_new);
        //const double dist = vector_distance(J_new, J);
        // cerr << "ITR: " << burnin_itr << ", J distance dev: " << dist << endl;
        //J = J_new;
        //D = D_new;
      }
      if (VERBOSE)
        cout << ", acc rate: " << n_acc / ((n_sites-2)*burnin) << endl;
      add_paths(paths, average_paths);
      write_root_to_pathfile_local(outfile, th.node_names.front());
      for (size_t node_id = 1; node_id < th.n_nodes; ++node_id)
        append_to_pathfile_local(outfile, th.node_names[node_id], paths[node_id]);
    }
    /* AVERAGING PATHS */
    if (n_samples > 0)
      for (size_t site_id = 0; site_id < average_paths.size(); site_id++)
        transform(average_paths[site_id].begin(), average_paths[site_id].end(),
                  average_paths[site_id].begin(),
                  std::bind(std::divides<double>(), std::placeholders::_1,
                            n_samples));
    
    /* (6) OUTPUT */
    if (VERBOSE)
      cerr << "[WRITING OUTPUT]" << endl;
    // write_output(outfile, average_paths);
        

  }

  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
