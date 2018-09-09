/* Copyright (C) 2018 University of Southern California
 *                    Jianghan Qu and Andrew D Smith
 *
 * Author: Andrew D. Smith and Jianghan Qu
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
#include <random>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <istream>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "PhyloTreePreorder.hpp"
#include "TreeHelper.hpp"
#include "TripletSampler.hpp"
#include "StateSeq.hpp"
#include "EpiEvoModel.hpp"
#include "GlobalJump.hpp"
#include "Path.hpp"
#include "ParamEstimation.hpp"

using std::vector;
using std::string;
using std::endl;
using std::cerr;
using std::cout;
using std::to_string;
using std::ostream_iterator;
using std::numeric_limits;
using std::istringstream;


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
    //iss.rdbuf()->pubsetbuf(const_cast<char*>(buffer.c_str()), buffer.size());
    iss.str(std::move(buffer));
    size_t site_index = 0;
    iss >> site_index; // not important info but must be removed

    // now read the states for the current site
    size_t node_idx = 0;
    T tmp_state_val;
    while (node_idx < n_nodes && iss >> tmp_state_val)
      state_sequences[node_idx++].push_back(tmp_state_val == '1');
    
    if (node_idx < n_nodes)
      throw std::runtime_error("inconsistent number of states: " +
                               to_string(node_idx) + "/" + to_string(n_nodes)
                               + "in line: " + to_string(site_index));

    ++site_count;
  }

  if (site_count == 0)
    throw std::runtime_error("no sites read from states file: " + statesfile);

  return site_count;
}


bool
file_is_readable(const string &param_file) {
  std::ifstream in(param_file.c_str());
  return in.good();
}


static void
write_output(const bool only_leaf_nodes, const TreeHelper &th,
             const string &outfile, const vector<string> &node_names,
             const vector<StateSeq> &sequences) {

  std::ofstream out(outfile.c_str());
  if (!out)
    throw std::runtime_error("bad output file: " + outfile);

  const size_t n_sequences = sequences.size();
  const size_t n_sites = sequences.front().seq.size();

  out << '#';
  bool first_name_written = false;
  for (size_t i = 0; i < th.n_nodes; ++i)
    if (!only_leaf_nodes || th.is_leaf(i)) {
      if (first_name_written)
        out << '\t';
      else first_name_written = true;
      out << node_names[i];
    }
  out << '\n';

  for (size_t i = 0; i < n_sites; ++i) {
    out << i;
    for (size_t j = 0; j < n_sequences; ++j)
      if (!only_leaf_nodes || th.is_leaf(j))
        out << '\t' << static_cast<bool>(sequences[j].seq[i]);
    out << '\n';
  }
}


/* This function does the sampling for an individual change in the
   state sequence
 */
static void
sample_jump(const EpiEvoModel &the_model, const double total_time,
            std::mt19937 &gen, TripletSampler &ts, vector<GlobalJump> &the_path,
            double &time_value) {

  static const size_t n_triplets = 8;

  // triplet_count = c_{ijk} for current sequence (encoded in the
  // TripletSampler object)
  vector<size_t> triplet_counts;
  ts.get_triplet_counts(triplet_counts);

  // holding_rate = c_{ijk}*lambda_{ijk}
  const double holding_rate =
    std::inner_product(triplet_counts.begin(), triplet_counts.end(),
                       the_model.triplet_rates.begin(), 0.0);

  // sample a holding time = time until next state change
  std::exponential_distribution<double> exp_distr(holding_rate);
  const double holding_time = std::max(exp_distr(gen),
                                       numeric_limits<double>::min());

  // update the current time_value
  time_value += holding_time;

  // if the holding time ends before the total time interval, we can
  // make a change to the state sequence
  if (time_value < total_time) {

    /* first: get a probability distribution for the triplet to change */
    vector<double> triplet_prob(n_triplets, 0.0);
    for (size_t i = 0; i < n_triplets; ++i)
      triplet_prob[i] =
        triplet_counts[i]*the_model.triplet_rates[i]/holding_rate;

    /* next: use that distribution to sample which triplet type at
       which the change will happen */
    std::discrete_distribution<size_t> multinom(triplet_prob.begin(),
                                                triplet_prob.end());
    const size_t context = multinom(gen);

    /* sample a change position having the relevant triplet; this
       changes the TripletSampler data structure to reflect a changed
       state at the position sampled */
    const size_t change_position = ts.random_mutate(context, gen);

    /* add the changed position and change time to the path */
    the_path.push_back(GlobalJump(time_value, change_position));
  }
}

static void
assign_reversed_changes_to_sites(const vector<GlobalJump> &global_path,
                        vector<Path> &by_site) {

  const size_t n_changes = global_path.size();
  for (size_t i = 0; i < n_changes; ++i) {
    const size_t j = n_changes - i - 1;
    const size_t pos = global_path[j].position;
    const double time = by_site[pos].tot_time - global_path[j].timepoint;
    by_site[pos].jumps.push_back(time);
  }
}


int main(int argc, const char **argv) {

  try {

    string outfile;
    string pathfile;
    string tree_file;
    string statfile;
    bool VERBOSE = false;
    bool write_only_leaves = false;
    size_t n_sites = 100;
    size_t n_samples = 1000;

    double evolutionary_time = numeric_limits<double>::lowest();

    size_t rng_seed = numeric_limits<size_t>::max();

    ////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]),
                           "collect sumamry statistics"
                           "from (stationary) simulation",
                           "<params-file>");
    opt_parse.add_opt("output", 'o', "name of output file for methylomes"
                      "(default: stdout)", false, outfile);
    opt_parse.add_opt("stat", 'S', "name of stat file", false, statfile);
    opt_parse.add_opt("n-sites", 'n', "length of sequence to simulate "
                      "(default: " + to_string(n_sites) + ")", false, n_sites);
    opt_parse.add_opt("n-samples", 'm', "number of histories to simulate "
                      "(default: " + to_string(n_samples) + ")",
                      false, n_samples);
    opt_parse.add_opt("paths", 'p', "name of output file for evolution paths"
                      "as sorted jump times (default: stdout)", false, pathfile);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);
    opt_parse.add_opt("tree", 't', "Newick format tree file", false, tree_file);
    opt_parse.add_opt("evo-time", 'T', "evolutionary time", false,
                      evolutionary_time);
    opt_parse.add_opt("leaf", 'l', "write only leaf states", false,
                      write_only_leaves);
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string param_file(leftover_args.front());
    if (!file_is_readable(param_file)) {
      cerr << "cannot read file: "<< param_file << endl;
      return  EXIT_SUCCESS;
    }
    if (!tree_file.empty()) {
      if (evolutionary_time != numeric_limits<double>::lowest()) {
        cerr << "specify exactly one of: tree or time" << endl;
        return  EXIT_SUCCESS;
      }
      if (!file_is_readable(tree_file)) {
        cerr << "cannot read file: "<< tree_file << endl;
        return  EXIT_SUCCESS;
      }
    }
    else if (evolutionary_time == numeric_limits<double>::lowest()) {
      cerr << "specify exactly one of: tree or time" << endl;
      return  EXIT_SUCCESS;
    }
    ////////////////////////////////////////////////////////////////////////

    /* (1) INITIALIZING PARAMETERS AND TREE */
    if (VERBOSE)
      cerr << "reading parameter file: " << param_file << endl;
    EpiEvoModel the_model;
    read_model(param_file, the_model);
    the_model.scale_triplet_rates();
    if (VERBOSE)
      cerr << the_model << endl;

    size_t n_nodes = 0;
    TreeHelper th;
    if (!tree_file.empty()) {
      if (VERBOSE)
        cerr << "reading tree file: " << tree_file << endl;
      PhyloTreePreorder the_tree; // tree topology and branch lengths
      std::ifstream tree_in(tree_file.c_str());
      if (!tree_in || !(tree_in >> the_tree))
        throw std::runtime_error("bad tree file: " + tree_file);
      n_nodes = the_tree.get_size();
      th = TreeHelper(the_tree);
    }
    else {
      if (VERBOSE)
        cerr << "initializing two node tree with time: "
             << evolutionary_time << endl;
      n_nodes = 2;
      th = TreeHelper(evolutionary_time);
    }

    /* standard mersenne_twister_engine seeded with rd()*/
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    if (VERBOSE)
      cerr << "rng seed: " << rng_seed << endl;
    std::mt19937 gen(rng_seed);

    /* (2) GET THE ROOT SEQUENCE */
    if (VERBOSE)
      cerr << "[OBTAINING ROOT SEQUENCE]" << endl;

    vector<char> root_seq;
    if (VERBOSE)
      cerr << "[SIMULATING: " << th.node_names[0] << " (ROOT)]" << endl;
    the_model.sample_state_sequence_init(n_sites, gen, root_seq);
    if (VERBOSE)
      cerr << "[ROOT LENGTH: " << n_sites << "]" << endl;
   
    StateSeq s(root_seq);
    if (VERBOSE) {
      cerr << "[SUMMARY:]" << endl
           << s.summary_string() << endl;
      vector<double> triplet_props;
      s.get_triplet_proportions(triplet_props);
      double total_rate = 0;
      for (size_t i = 0; i < triplet_props.size(); ++i)
        total_rate += the_model.triplet_rates[i]*triplet_props[i];
      // do we need to divide by the sequene length here?
      cerr << "mutations per site (at root): " << total_rate << endl;
    }

    if (!pathfile.empty())
      write_root_to_pathfile_global(pathfile, th.node_names[0], s);

    vector<StateSeq> sequences(n_nodes, s);

    /* HEADER LINE OF OUTPUT STAT FILE */
    std::ofstream of;
    of.open(statfile.c_str());
    std::ostream out(of.rdbuf());
    if (!out)
      throw std::runtime_error("bad output file: " + statfile);
    out << "J_000\tJ_001\tJ_010\tJ_011\tJ_100\tJ_101\tJ_110\tJ_111\t"
    << "D_000\tD_001\tD_010\tD_011\tD_100\tD_101\tD_110\tD_111\n" << endl;
  
    
    for (size_t sample_id = 0; sample_id < n_samples; sample_id++) {
      /* SIMULATION */
      vector<vector<Path> > all_paths;
      
      // root node
      all_paths.push_back(vector<Path> (n_sites, Path()));
      
      // iterate over branches
      for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
        vector<Path> branch_paths;
        const double curr_branch_len = th.branches[node_id];
        
        TripletSampler ts(sequences[th.parent_ids[node_id]]);
        double time_value = 0;
        vector<GlobalJump> the_path;
        
        while (time_value < curr_branch_len)
          sample_jump(the_model, curr_branch_len, gen, ts, the_path, time_value);
        
        ts.get_sequence(sequences[node_id]);
        
        // Convert REVERTED global jumps to local paths
        vector<Path> path_by_site(n_sites);
        for (size_t j = 0; j < path_by_site.size(); ++j) {
          path_by_site[j].init_state = sequences[node_id].seq[j];
          path_by_site[j].tot_time = curr_branch_len;
        }
        assign_reversed_changes_to_sites(the_path, path_by_site);
        // Now we sepcifically want to test branch 1 only
        if (node_id == 1)
          all_paths.push_back(path_by_site);
                        
        if (sample_id == 0 && !pathfile.empty()) {
          // reverse global jumps
          vector<GlobalJump> reversed_jumps(the_path);
          std::reverse(reversed_jumps.begin(), reversed_jumps.end());
          for(size_t j = 0; j < reversed_jumps.size(); j++) {
            reversed_jumps[j].timepoint = curr_branch_len -
            reversed_jumps[j].timepoint;
          }

          // only output one instance of simualted paths
          append_to_pathfile_global(pathfile, th.node_names[node_id],
                                    reversed_jumps);
          if (!outfile.empty())
            write_output(write_only_leaves, th, outfile, th.node_names,
                         sequences);
        }
      }
      
      /* EXTRACT SUMMARY STATISTICS */
      vector<double> J;
      vector<double> D;
      get_sufficient_statistics(all_paths, J, D);
      for (size_t i = 0; i < 8; i++)
        out << J[i] << "\t";
      for (size_t i = 0; i < 7; i++)
        out << D[i] << "\t";
      out << D[7] << endl;

    }

  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
