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
#include <array>
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
using std::array;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::istringstream;
using std::ofstream;
using std::runtime_error;
using std::to_string;
using std::numeric_limits;
using std::ostream_iterator;
using std::ifstream;


template <class T>
size_t
read_states_file(const string &statesfile, vector<vector<T> > &state_sequences,
                 const TreeHelper &th) {

  std::ifstream in(statesfile);
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
    else if (th.is_leaf(node_id))
      throw std::runtime_error("no data in leaf node: " + node_name);
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
template <typename RandEngine,
          typename StateType>
static void
initialize_paths(RandEngine &gen, const TreeHelper &th,
                 vector<vector<StateType> > &state_sequences,
                 vector<vector<Path> > &paths) {

  const size_t n_sites = state_sequences.front().size();

  paths.resize(n_sites);
  for (size_t site_id = 0; site_id < n_sites; ++site_id) {
    paths[site_id].clear();
    paths[site_id].resize(th.n_nodes);
  }

  auto unif =
    bind(std::uniform_real_distribution<double>(0.0, 1.0), std::ref(gen));

  vector<StateType> child_states = {0, 0}; // two child states

  for (size_t i = th.n_nodes; i > 0; --i) {

    const size_t node_id = i - 1;

    for (size_t site_id = 0; site_id < n_sites; ++site_id) {

      if (!th.is_leaf(node_id)) {
        // get child states
        size_t n_ch = 0;
        for (ChildSet c(th.subtree_sizes, node_id); c.good(); ++c)
          child_states[n_ch++] = state_sequences[*c][site_id];

        // sample parent state
        if (th.is_root(node_id))
          state_sequences[node_id][site_id] = state_sequences[0][site_id];
        else
          state_sequences[node_id][site_id] = child_states[std::floor(unif()*n_ch)];

        // assign site-specific paths above two child nodes
        for (ChildSet c(th.subtree_sizes, node_id); c.good(); ++c) {

          const double len = th.branches[*c];
          paths[site_id][*c].tot_time = len;
          paths[site_id][*c].init_state = state_sequences[node_id][site_id];

          if (state_sequences[*c][site_id] != paths[site_id][*c].init_state)
            paths[site_id][*c].jumps.push_back(unif() * len);
        }
      }
    }
  }
}


static void
sample_summary_stats(const vector<double> &rates, const TreeHelper &th,
                     vector<vector<Path> > &paths,
                     vector<vector<double> > &J_all_sites,
                     vector<vector<double> > &D_all_sites,
                     std::mt19937 &gen,
                     const size_t batch) {
  const size_t n_triples = 8;
  J_all_sites.resize(th.n_nodes);
  D_all_sites.resize(th.n_nodes);
  for(size_t b = 1; b < th.n_nodes; ++b) {
    J_all_sites[b].clear();
    J_all_sites[b].resize(n_triples, 0.0);
    D_all_sites[b].clear();
    D_all_sites[b].resize(n_triples, 0.0);
  }

  for (size_t i = 0; i < batch; i++) {
    update_paths_indep(rates, th, paths, gen);

    vector<vector<double> > J_one_site, D_one_site;
    get_sufficient_statistics(paths, J_one_site, D_one_site);
    for (size_t b = 1; b < th.n_nodes; b++) {
      for (size_t i = 0; i < n_triples; i++) {
        J_all_sites[b][i] += J_one_site[b][i];
        D_all_sites[b][i] += D_one_site[b][i];
      }
    }
  }

  /* CALCULATE BATCH AVERAGE */
  for (size_t b = 1; b < th.n_nodes; ++b) {
    for (size_t i = 0; i < n_triples; i++) {
      J_all_sites[b][i] /= batch;
      D_all_sites[b][i] /= batch;
    }
  }
}



static void
initialize_model_from_indep_rates(EpiEvoModel &the_model,
                                  const vector<double> rates) {

  the_model.stationary_baseline.reset();
  the_model.T.reset();
  the_model.Q.reset();

  // set r_0_ = r0, r_1_ = r1
  array<double, 8> triplet_rates;
  for (size_t i = 0; i < the_model.n_triplets; i++)
    triplet_rates[i] = rates[(i/2) % 2];

  the_model.rebuild_from_triplet_rates(triplet_rates);
}


////////////////////////////////////////////////////////////////////////////////
///////////////                  output                /////////////////////////
////////////////////////////////////////////////////////////////////////////////

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

static void
write_mcmc_verbose_header(std::ostream &out) {
  vector<string> header_tokens = {
                                  "itr",
                                  "rate0",
                                  "rate1",
                                  "tree",
  };
  copy(begin(header_tokens), end(header_tokens),
       std::ostream_iterator<string>(out, "\t"));
  out << '\n';
}


int main(int argc, const char **argv) {

  try {

    static const double param_tol = 1e-10;

    bool VERBOSE = false;
    bool optimize_branches = false;
    double evolutionary_time = 0.0;

    size_t rng_seed = numeric_limits<size_t>::max();
    size_t iterations = 10;
    size_t batch = 10;              // MCMC iterations

    string paramfile;
    string pathfile;
    string tree_file, treefile_updated;
    ////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]),
                           "generate initial paths and parameters"
                           " given states at leaves",
                           "(<tree-file>) <states-file>");
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);
    opt_parse.add_opt("iterations", 'i', "number of iterations",
                      false, iterations);
    opt_parse.add_opt("batch", 'B', "number of MCMC iteration", false, batch);
    opt_parse.add_opt("param", 'p', "output file of parameters",
                      false, paramfile);
    opt_parse.add_opt("outtree", 't', "output file of tree",
                      false, treefile_updated);
    opt_parse.add_opt("evo-time", 'T', "evolutionary time (assumes no tree)",
                      false, evolutionary_time);
    opt_parse.add_opt("path", 'o', "output file of local paths (default: stdout)",
                      false, pathfile);
    opt_parse.add_opt("branch", 'b', "optimize branch lengths as well",
                      false, optimize_branches);
    opt_parse.set_show_defaults();
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
    if (leftover_args.size() == 1) {
      if (evolutionary_time == 0.0) {
        cerr << opt_parse.help_message() << endl;
        return EXIT_SUCCESS;
      }
    }
    else if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    else {
      tree_file = leftover_args.front();
    }
    const string statesfile(leftover_args.back());
    ////////////////////////////////////////////////////////////////////////

    PhyloTreePreorder the_tree; // tree topology and branch lengths
    TreeHelper th;
    if (evolutionary_time > 0.0) {
      if (VERBOSE)
        cerr << "[INITIALIZING TWO NODE TREE WITH TIME: "
             << evolutionary_time << "]" << endl;
      th = TreeHelper(evolutionary_time);
    }
    else {
      if (VERBOSE)
        cerr << "[READING TREE: " << tree_file << "]" << endl;
      ifstream tree_in(tree_file);
      if (!tree_in || !(tree_in >> the_tree))
        throw runtime_error("bad tree file: " + tree_file);
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
    vector<vector<Path> > paths; // [sites] x [nodes]
    initialize_paths(gen, th, state_sequences, paths);

    /* Run EM to learn a site-independent model */
    vector<double> rates(2, 0.0);
    vector<vector<double> > J, D;
    compute_sufficient_statistics(paths, J, D);

    if (VERBOSE) {
      write_mcmc_verbose_header(cerr);
      cerr << "0" << "\t" << rates[0] << "\t" << rates[1] << endl;
    }

    for (size_t itr = 0; itr < iterations; itr++) {
      if (!optimize_branches)
        estimate_rates_indep(J, D, rates, th);
      else {
        estimate_rates_and_branches_indep(J, D, rates, th, paths);
        the_tree.set_branch_lengths(th.branches);
      }

      expectation_sufficient_statistics(rates, th, paths, J, D);
      // Report
      if (VERBOSE)
        cerr << itr+1 << "\t" << rates[0] << "\t" << rates[1] << endl;
    }

    /* Re-sample a better initial path */
    vector<vector<double> > J_trip, D_trip;
    sample_summary_stats(rates, th, paths, J_trip, D_trip, gen, batch);

    /*******************************************************/
    /* Generate initial parameters of context-dependent model */
    /*******************************************************/

    if (VERBOSE)
      cerr << "[CONSTRUCTING EPIEVO MODEL]" << endl;

    // initialize a EpievoModel
    EpiEvoModel the_model;
    initialize_model_from_indep_rates(the_model, rates);

    // re-estimate triplet rates from paths
    if (!optimize_branches) {
      estimate_rates(false, param_tol, J_trip, D_trip, the_model);
      set_one_change_per_site_per_unit_time(the_model.triplet_rates,
                                            th.branches);
    }
    else {
      estimate_rates_and_branches(false, param_tol, J_trip, D_trip, th,
                                  the_model);
      the_tree.set_branch_lengths(th.branches);
    }
    scale_jump_times(paths, th);

    // write path file
    if (VERBOSE)
      cerr << "[WRITING PATHS]" << endl;

    write_root_to_pathfile_local(pathfile, th.node_names.front());
    for (size_t node_id = 1; node_id < th.n_nodes; ++node_id)
      append_to_pathfile_local(pathfile, th.node_names[node_id],
                               paths, node_id);

    // write parameters
    if (VERBOSE)
      cerr << "[WRITING PARAMETERS]\n" << the_model << endl;
    if (VERBOSE && evolutionary_time == 0.0)
      cerr << the_tree << endl;

    std::ofstream of_param;
    if (!paramfile.empty()) of_param.open(paramfile.c_str());
    std::ostream out_param(paramfile.empty() ?
                           std::cout.rdbuf() : of_param.rdbuf());
    if (!out_param)
      throw std::runtime_error("bad output param file: " + paramfile);

    out_param << the_model.format_for_param_file() << endl;

    if (optimize_branches) {
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
