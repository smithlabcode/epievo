/* Copyright (C) 2019 University of Southern California
 *                    Xiaojing Ji, and Andrew D Smith
 *
 * Author: Andrew D. Smith and Xiaojing Ji
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
 */

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <functional>
#include <algorithm>
#include <sstream>
#include <exception>
#include <random>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "PhyloTreePreorder.hpp"
#include "Path.hpp"
#include "TreeHelper.hpp"
#include "IndepSite.hpp"
#include "ParamEstimation.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::istringstream;
using std::runtime_error;
using std::to_string;
using std::numeric_limits;
using std::ostream_iterator;



template <class T>
size_t
read_states_file(const string &statesfile, vector<vector<T> > &state_sequences,
                 const TreeHelper &th) {
  
  std::ifstream in(statesfile.c_str());
  if (!in)
    throw std::runtime_error("bad states file: " + statesfile);

  // first line is the list of node names
  string buffer;
  if (!getline(in, buffer))
    throw std::runtime_error("cannot read nodes line in: " + statesfile);

  std::istringstream nodes_iss(buffer);
  vector<string> node_names_in;
  
  string tmp_node_name;
  while (nodes_iss >> tmp_node_name)
    node_names_in.push_back(tmp_node_name);

  if (node_names_in.size() < 2)
    throw std::runtime_error("fewer than 2 nodes names in: " + statesfile);

  if (node_names_in.front()[0] == '#') {
    if (node_names_in.front().length() == 1)
      node_names_in = vector<string>(node_names_in.begin() + 1,
                                      node_names_in.end());
    else node_names_in.front() = node_names_in.front().substr(1);
  }
  
  // get sorting indices of nodes in pre-order
  vector<size_t> idx_in_tree (node_names_in.size(), th.n_nodes);

  for (size_t node_id = 0; node_id < th.n_nodes; node_id++) {
    const string node_name = th.node_names[node_id];
    vector<string>::iterator it = std::find(node_names_in.begin(),
                                            node_names_in.end(),
                                            node_name);
    if (it != node_names_in.end())
      idx_in_tree[std::distance(node_names_in.begin(), it)] = node_id;
    else
      throw std::runtime_error("no data in node: " + node_name);
  }

  // read states
  const size_t n_nodes_in = node_names_in.size();
  size_t site_count = 0;

  state_sequences.resize(th.n_nodes);
  while (getline(in, buffer)) {
    istringstream iss;
    iss.str(std::move(buffer));

    size_t site_index = 0;
    iss >> site_index; // not important info but must be removed
    // now read the states for the current site
    size_t node_idx = 0;
    T tmp_state_val;
    while (node_idx < n_nodes_in && iss >> tmp_state_val)
      if (idx_in_tree[node_idx] < th.n_nodes)
        state_sequences[idx_in_tree[node_idx++]].push_back(tmp_state_val);
    if (node_idx < n_nodes_in)
      throw std::runtime_error("inconsistent number of states: " +
                               to_string(node_idx) + "/" +
                               to_string(n_nodes_in));

    ++site_count;
  }

  if (site_count == 0)
    throw std::runtime_error("no sites read from states file: " + statesfile);

  // resize internal state sequences
  for (size_t node_id = 0; node_id < th.n_nodes; node_id++)
    if (state_sequences[node_id].size() < site_count)
      state_sequences[node_id].resize(site_count);
  
  return site_count;
}


////////////////////////////////////////////////////////////////////////////////
///////////////              initialization            /////////////////////////
////////////////////////////////////////////////////////////////////////////////

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
        if (!th.is_root(node_id))
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
initialize_model_from_indep_rates(EpiEvoModel &the_model,
                                  const vector<double> rates) {
  const vector<vector<double> > zeros_two_by_two =
  vector<vector<double> > (2, vector<double> (2, 0.0));
  
  the_model.stationary_logbaseline = zeros_two_by_two;
  the_model.T = zeros_two_by_two;
  the_model.init_T = zeros_two_by_two;
  the_model.Q = zeros_two_by_two;
  
  // set r_0_ = r0, r_1_ = r1
  vector<double> triplet_rates (the_model.n_triplets);
  for (size_t i = 0; i < the_model.n_triplets; i++)
    triplet_rates[i] = rates[i / 2 % 2];
  
  the_model.rebuild_from_triplet_rates(triplet_rates);
}


////////////////////////////////////////////////////////////////////////////////
///////////////               write output             /////////////////////////
////////////////////////////////////////////////////////////////////////////////

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


int main(int argc, const char **argv) {

  try {

    static const double param_tol = 1e-10;
    
    bool VERBOSE = false;
    bool OPTBRANCH = false;
    bool ONEBRANCH = false;

    size_t rng_seed = numeric_limits<size_t>::max();
    size_t iterations = 10;
    
    string paramfile;
    string pathfile;
    string treefile_updated;
    
    ////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]),
                           "generate initial paths and parameters"
                           " given states at leaves",
                           " <tree> <states>");
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);
    opt_parse.add_opt("iterations", 'i', "number of iterations",
                      false, iterations);
    opt_parse.add_opt("param", 'p',
                      "output file of parameters (default: stdout)",
                      false, paramfile);
    opt_parse.add_opt("outtree", 't',
                      "output file of tree (default: stdout)",
                      false, treefile_updated);
    opt_parse.add_opt("one-branch", 'T', "one-branch tree", false,
                      ONEBRANCH);
    opt_parse.add_opt("path", 'o', "output file of local paths (default: stdout)",
                      false, pathfile);
    opt_parse.add_opt("branch", 'b', "optimize branch lengths as well",
                      false, OPTBRANCH);

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
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string tree(leftover_args[0]);
    const string statesfile(leftover_args[1]);
    ////////////////////////////////////////////////////////////////////////

    if (VERBOSE)
      cerr << "[READING TREE: " << tree << "]" << endl;
    PhyloTreePreorder the_tree; // tree topology and branch lengths
    TreeHelper th;
    if (ONEBRANCH) {
      if (VERBOSE)
        cerr << "initializing two node tree with time: " << tree << endl;
      th = TreeHelper(std::stod(tree));
    } else {
      cerr << "reading tree file: " << tree << endl;
      std::ifstream tree_in(tree.c_str());
      if (!tree_in || !(tree_in >> the_tree))
        throw std::runtime_error("bad tree file: " + tree);
      th = TreeHelper(the_tree);
    }

    vector<vector<bool> > state_sequences;

    if (VERBOSE)
      cerr << "[READING STATES FILE: " << statesfile << "]" << endl;
    read_states_file(statesfile, state_sequences, th);
    
    
    /* standard mersenne_twister_engine seeded with rd()*/
    if (rng_seed == numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    if (VERBOSE)
      cerr << "rng seed: " << rng_seed << endl;
    std::mt19937 gen(rng_seed);
    
    /* generate initial paths by heuristics */
    vector<vector<Path> > paths; // along multiple branches
    initialize_paths(gen, th, state_sequences, paths);

    /* Run EM to learn a site-independent model */
    vector<double> rates (2, 0.0);
    vector<double> init_pi (2, 0.0);
    vector<double> init_pi_post (2, 0.0);
    vector<vector<double> > J, D;
    compute_sufficient_statistics(paths, J, D);
    estimate_root_distribution(paths, init_pi);
 
    if (VERBOSE)
      cerr << "ITR\tRATE0\tRATE1\t\tINIT0\t\tINIT1\tTREE\n"
      << "0" << "\t" << rates[0] << "\t" << rates[1] << "\t"
      << init_pi[0] << "\t" << init_pi[1] << endl;

    
    for (size_t itr = 0; itr < iterations; itr++) {
      if (!OPTBRANCH)
        estimate_rates(J, D, rates, th);
      else {
        estimate_rates_and_branches(J, D, rates, th, paths);
        the_tree.set_branch_lengths(th.branches);
      }

      expectation_sufficient_statistics(rates, init_pi, th, paths, J, D,
                                        init_pi_post);

      init_pi = init_pi_post;
      
      // Report
      if (VERBOSE)
        cerr << itr+1 << "\t" << rates[0] << "\t" << rates[1] << "\t"
        << init_pi[0] << "\t" << init_pi[1] << endl;
    }
    
    /* Re-sample a better initial path */
    vector<vector<Path> > sampled_paths(paths);
    sample_paths(rates, init_pi, th, paths, gen, sampled_paths);
    
    /*******************************************************/
    /* Generate initial parameters of context-dependent model */
    /*******************************************************/

    if (VERBOSE)
      cerr << "[CONSTRUCTING EPIEVO MODEL]" << endl;

    // initialize a EpievoModel
    EpiEvoModel the_model;
    initialize_model_from_indep_rates(the_model, rates);

    // set init_T
    estimate_root_distribution(sampled_paths, the_model);
    
    // re-estimate triplet rates from paths
    if (!OPTBRANCH)
      compute_estimates_for_rates_only(false, param_tol, sampled_paths,
                                       the_model);
    else {
      compute_estimates_rates_and_branches(false, param_tol, sampled_paths,
                                           th, the_model);
      scale_jump_times(sampled_paths, th);
      the_tree.set_branch_lengths(th.branches);
    }

    // write path file
    if (VERBOSE)
      cerr << "[WRITING PATHS]" << endl;

    write_root_to_pathfile_local(pathfile, th.node_names.front());
    for (size_t node_id = 1; node_id < th.n_nodes; ++node_id)
      append_to_pathfile_local(pathfile, th.node_names[node_id],
                               sampled_paths[node_id]);
    
    // write parameters
    if (VERBOSE)
      cerr << "[WRITING PARAMETERS]\n" << the_model << "\n" << the_tree << endl;

    std::ofstream of_param;
    if (!paramfile.empty()) of_param.open(paramfile.c_str());
    std::ostream out_param(paramfile.empty() ?
                           std::cout.rdbuf() : of_param.rdbuf());
    if (!out_param)
      throw std::runtime_error("bad output param file: " + paramfile);
    
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
