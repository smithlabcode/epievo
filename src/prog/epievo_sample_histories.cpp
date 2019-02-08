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
#include <sstream>
#include <exception>
#include <random>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "ParamEstimation.hpp"
#include "PhyloTreePreorder.hpp"
#include "Path.hpp"
#include "StateSeq.hpp"
#include "EpiEvoModel.hpp"
#include "TreeHelper.hpp"
#include "SingleSiteSampler.hpp"

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

static void
write_root_to_pathfile_local(const string &outfile, const string &root_name) {
  std::ofstream outpath(outfile.c_str());
  outpath << "NODE:" << root_name << endl;
}


template <class T>
static void
write_output(const string &outfile, const vector<string> &node_names,
             const vector<vector<T> > &sequences) {

  std::ofstream out(outfile.c_str());
  if (!out)
    throw std::runtime_error("bad output file: " + outfile);

  const size_t n_sites = sequences.front().size();
  const size_t n_nodes = sequences.size();

  out << '#';
  copy(begin(node_names), end(node_names), ostream_iterator<string>(out, "\t"));
  out << '\n';

  for (size_t i = 0; i < n_sites; ++i) {
    out << i;
    for (size_t j = 0; j < n_nodes; ++j)
      out << '\t' << sequences[j][i];
    out << '\n';
  }
}



template <class T>
size_t
read_states_file(const string &statesfile,
                 vector<string> &node_names,
                 vector<vector<T> > &state_sequences) {

  std::ifstream in(statesfile.c_str());
  if (!in)
    throw std::runtime_error("bad states file: " + statesfile);

  // first line should give the node names in pre-order traversal
  string buffer;
  if (!getline(in, buffer))
    throw std::runtime_error("cannot read nodes line in: " + statesfile);

  std::istringstream nodes_iss(buffer);
  node_names.clear();
  string tmp_node_name;
  while (nodes_iss >> tmp_node_name)
    node_names.push_back(tmp_node_name);

  if (node_names.size() < 2)
    throw std::runtime_error("fewer than 2 nodes names in: " + statesfile);

  // remove leading '#' from first node name if it exists
  if (node_names.front()[0] == '#') {
    if (node_names.front().length() == 1)
      node_names = vector<string>(node_names.begin() + 1, node_names.end());
    else node_names.front() = node_names.front().substr(1);
  }
  const size_t n_nodes = node_names.size();
  state_sequences = vector<vector<T> >(node_names.size());

  size_t site_count = 0;
  while (getline(in, buffer)) {
    istringstream iss;
    iss.rdbuf()->pubsetbuf(const_cast<char*>(buffer.c_str()), buffer.size());

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


template <class T>
static void
initialize_internal_states(std::mt19937 &gen,
                           const TreeHelper &th,
                           vector<vector<T> > &state_sequences) {

  const size_t n_sites = state_sequences.front().size();

  auto unif =
    bind(std::uniform_real_distribution<double>(0.0, 1.0), std::ref(gen));

  vector<T> child_states = {0, 0};
  for (size_t site_id = 0; site_id < n_sites; ++site_id) {
    for (size_t node_id = th.n_nodes; node_id > 0; --node_id) {
      cerr << site_id << '\t' << node_id - 1 << endl;
      if (!th.is_leaf(node_id-1)) {
        size_t n_ch = 0;
        for (auto j = ChildSet(th.subtree_sizes, node_id-1); j.good(); ++j)
          child_states[n_ch++] = state_sequences[*j][site_id];
        state_sequences[node_id-1][site_id] =
          child_states[std::floor(unif()*n_ch)];
      }
    }
  }
}


static void
append_to_pathfile_local(const string &pathfile, const string &node_name,
                         const vector<Path> &path_by_site) {
  std::ofstream outpath(pathfile.c_str(), std::ofstream::app);
  outpath << "NODE:" << node_name << endl;
  for (size_t i = 0; i < path_by_site.size(); ++i)
    outpath << i << '\t' << path_by_site[i] << '\n';
}


int main(int argc, const char **argv) {

  try {

    static const double param_tol = 1e-10;

    string outfile;
    bool VERBOSE = false;
    bool input_includes_internal = false;
    size_t rng_seed = numeric_limits<size_t>::max();
    size_t n_samples = numeric_limits<size_t>::max();
    bool input_is_paths = false;

    ////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]),
                           "sample histories given starting data"
                           " as histories or states (at leaves or internal)",
                           "<params> <tree> <states/paths>");
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);
    opt_parse.add_opt("internal", 'i', "input includes internal states",
                      false, input_includes_internal);
    opt_parse.add_opt("n-samples", 'n', "number of samples",
                      false, n_samples);
    opt_parse.add_opt("paths", 'P', "input is paths",
                      false, input_is_paths);
    opt_parse.add_opt("output", 'o', "output parameter file",
                      false, outfile);
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
    if (leftover_args.size() != 3) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string param_file(leftover_args[0]);
    const string treefile(leftover_args[1]);
    const string input_file(leftover_args[2]);
    ////////////////////////////////////////////////////////////////////////

    if (VERBOSE)
      cerr << "[READING TREE: " << treefile << "]" << endl;
    PhyloTreePreorder the_tree; // tree topology and branch lengths
    std::ifstream tree_in(treefile.c_str());
    if (!tree_in || !(tree_in >> the_tree))
      throw runtime_error("cannot read tree file: " + treefile);
    const TreeHelper th(the_tree);

    if (VERBOSE)
      cerr << "[READING PARAMETERS: " << param_file << endl;
    EpiEvoModel the_model;
    read_model(param_file, the_model);
    the_model.scale_triplet_rates();
    if (VERBOSE)
      cerr << the_model << endl;

    vector<string> node_names;
    vector<vector<Path> > paths; // along multiple branches
    vector<vector<bool> > state_sequences;
    size_t n_sites = 0;
    if (input_is_paths) {
      if (VERBOSE)
        cerr << "[READING PATHS FILE: " << input_file << "]" << endl;
      read_paths(input_file, node_names, paths);
      n_sites = paths.back().size();
    }
    else {
      if (VERBOSE)
        cerr << "[READING STATES FILE: " << input_file
             << (input_includes_internal ?
                 " (WITH INTERNAL)]\n" : " (LEAF ONLY)]\n");
      n_sites = read_states_file(input_file, node_names,
                                 state_sequences);
    }

    if (VERBOSE)
      cerr << "[CHECKING CONSISTENCY OF TREE AND STATES]" << endl
           << "node names: tree / states" << endl;
    size_t states_names_idx = 0;
    for (size_t i = 0; i < th.n_nodes; ++i)
      if (input_includes_internal || th.is_leaf(i)) {
        if (VERBOSE)
          cerr << th.node_names[i] << " / "
               << node_names[states_names_idx] << endl;
        if (th.node_names[i] != node_names[states_names_idx++])
          throw runtime_error("inconsistent node names in states vs tree file");
      }

    /* standard mersenne_twister_engine seeded with rd()*/
    if (rng_seed == numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    if (VERBOSE)
      cerr << "rng seed: " << rng_seed << endl;
    std::mt19937 gen(rng_seed);

    if (!input_is_paths) {
      if (!input_includes_internal) {
        vector<vector<bool> > tmp_state_seqs;
        tmp_state_seqs.swap(state_sequences);
        state_sequences.resize(th.n_nodes);
        size_t leaf_idx = 0;
        for (size_t i = 0; i < th.n_nodes; ++i) {
          if (th.is_leaf(i))
            swap(state_sequences[i], tmp_state_seqs[leaf_idx++]);
          else state_sequences[i].resize(n_sites);
        }
        initialize_internal_states(gen, th, state_sequences);
      }
      // nodes X sites
      paths = vector<vector<Path> >(th.n_nodes, vector<Path>(n_sites));
      for (size_t i = 1; i < th.n_nodes; ++i)
        for (size_t j = 0; j < n_sites; ++j)
          paths[i][j] =
            Path(state_sequences[th.parent_ids[i]][j], th.branches[i]);
    }

    if (n_samples == numeric_limits<size_t>::max())
      n_samples = n_sites;

    auto site_sampler =
      bind(std::uniform_int_distribution<size_t>(1, n_sites-2), std::ref(gen));

    for (size_t i = 0; i < n_samples; ++i) {
      size_t site_id = site_sampler();
      vector<Path> sampled_path;
      Metropolis_Hastings_site(the_model, th, site_id, paths, gen, sampled_path);
    }

    write_root_to_pathfile_local(outfile, th.node_names.front());
    for (size_t node_id = 1; node_id < th.n_nodes; ++node_id)
      append_to_pathfile_local(outfile, th.node_names[node_id], paths[node_id]);
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
