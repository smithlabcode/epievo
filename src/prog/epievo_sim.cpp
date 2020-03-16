/* Copyright (C) 2019 University of Southern California
 *                    Jianghan Qu, Andrew D Smith and Xiaojing Ji
 *
 * Author: Andrew D. Smith, Xiaojing Ji and Jianghan Qu
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
#include <bitset>


#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "PhyloTreePreorder.hpp"
#include "TreeHelper.hpp"
#include "TripletSampler.hpp"
#include "EpiEvoModel.hpp"
#include "GlobalJump.hpp"

#include "epievo_utils.hpp"

using std::vector;
using std::string;
using std::endl;
using std::cerr;
using std::cout;
using std::to_string;
using std::ostream_iterator;
using std::numeric_limits;
using std::istringstream;
using std::ifstream;


bool
file_is_readable(const string &param_file) {
  std::ifstream in(param_file);
  return in.good();
}


static void
write_output(const bool only_leaf_nodes, const TreeHelper &th,
             const string &outfile, const vector<string> &node_names,
             const vector<state_seq> &sequences) {

  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("bad output file: " + outfile);

  const size_t n_sequences = sequences.size();
  const size_t n_sites = sequences.front().size();

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
        out << '\t' << static_cast<bool>(sequences[j][i]);
    out << '\n';
  }
}


/* This function does the sampling for an individual change in the
   state sequence
 */
static void
sample_jump(const EpiEvoModel &the_model, const double total_time,
            std::mt19937 &gen, TripletSampler &ts, vector<GlobalJump> &the_path,
            double &time_value, vector<size_t> &events) {

  static const size_t n_triplets = 8;

  // triplet_count = c_{ijk} for current sequence (encoded in the
  // TripletSampler object)
  vector<size_t> triplet_counts;
  ts.get_triplet_counts(triplet_counts);

  // holding_rate = c_{ijk}*lambda_{ijk}
  const double holding_rate =
    inner_product(begin(triplet_counts), end(triplet_counts),
                  begin(the_model.triplet_rates), 0.0);

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
    std::discrete_distribution<size_t> multinom(begin(triplet_prob),
                                                end(triplet_prob));
    const size_t context = multinom(gen);
    ++events[context];

    /* sample a change position having the relevant triplet; this
       changes the TripletSampler data structure to reflect a changed
       state at the position sampled */
    const size_t change_position = ts.random_mutate(context, gen);

    /* add the changed position and change time to the path */
    the_path.push_back(GlobalJump(time_value, change_position));
  }
}


int main(int argc, const char **argv) {

  try {

    string pathfile;
    string tree_file;
    string root_states_file;
    bool VERBOSE = false;
    bool unscaled_model_params = false;
    bool scale_time = false;
    bool TRPARAM = false;
    bool write_only_leaves = false;
    size_t n_sites = 100;

    double evolutionary_time = numeric_limits<double>::lowest();

    size_t rng_seed = numeric_limits<size_t>::max();

    ////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]), "simulate epigenome evolution",
                           "<params-file> <outfile>");
    opt_parse.add_opt("n-sites", 'n', "length of sequence to simulate",
                      false, n_sites);
    opt_parse.add_opt("paths", 'p', "name of output file for evolution paths "
                      "as sorted jump times", false, pathfile);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);
    opt_parse.add_opt("root", 'r', "root states file", false, root_states_file);
    opt_parse.add_opt("tree", 't', "Newick format tree file", false, tree_file);
    opt_parse.add_opt("evo-time", 'T', "evolutionary time", false,
                      evolutionary_time);
    opt_parse.add_opt("leaf", 'l', "write only leaf states (default: all nodes)",
                      false, write_only_leaves);
    opt_parse.add_opt("unscaled-param", '\0', "do not scale model parameters",
                      false, unscaled_model_params);
    opt_parse.add_opt("scale-time", '\0', "scale time", false, scale_time);
    opt_parse.add_opt("rates", 'R', "use triplet transition rates",
                      false, TRPARAM);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
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
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string param_file(leftover_args.front());
    const string outfile(leftover_args.back());
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
    if (scale_time) {
      const double factor = rate_scaling_factor(the_model.triplet_rates);
      evolutionary_time /= factor;
    }
    if (!unscaled_model_params)
      the_model.scale_triplet_rates();

    if (VERBOSE)
      cerr << the_model << endl
           << the_model.format_for_param_file() << endl;

    size_t n_nodes = 0;
    TreeHelper th;
    if (!tree_file.empty()) {
      if (VERBOSE)
        cerr << "reading tree file: " << tree_file << endl;
      PhyloTreePreorder the_tree; // tree topology and branch lengths
      ifstream tree_in(tree_file);
      if (!tree_in || !(tree_in >> the_tree))
        throw std::runtime_error("bad tree file: " + tree_file);
      n_nodes = the_tree.get_size();
      th = TreeHelper(the_tree);
    }
    else {
      if (VERBOSE)
        cerr << "[initializing two node tree with time: "
             << evolutionary_time << "]" << endl;
      n_nodes = 2;
      th = TreeHelper(evolutionary_time);
    }

    /* standard mersenne_twister_engine seeded with rd()*/
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    if (VERBOSE)
      cerr << "[rng seed: " << rng_seed << "]" << endl;
    std::mt19937 gen(rng_seed);

    // get the root sequence
    if (VERBOSE)
      cerr << "[OBTAINING ROOT SEQUENCE]" << endl;

    state_seq root_seq;
    if (root_states_file.empty()) {
      if (VERBOSE)
        cerr << "[SIMULATING: " << th.node_names[0] << " (ROOT)]" << endl;
      the_model.sample_state_sequence_init(n_sites, gen, root_seq);
      if (VERBOSE)
        cerr << "[ROOT LENGTH: " << n_sites << "]" << endl;
    }
    else {
      if (VERBOSE)
        cerr << "[READING ROOT FILE: " << root_states_file << "]" << endl;
      vector<string> node_names;
      vector<state_seq> state_sequences;
      read_states_file(root_states_file, node_names, state_sequences);
      root_seq = state_sequences.front();
      n_sites = root_seq.size();
      if (VERBOSE)
        cerr << "[LOADED: " << th.node_names[0] << " (ROOT)]" << endl
             << "[ROOT LENGTH: " << n_sites << "]" << endl;
    }

    state_seq s(root_seq);
    if (VERBOSE) {
      cerr << "[SUMMARY:]" << endl
           << summary_string(s) << endl;
      vector<double> triplet_props;
      get_triplet_proportions(s, triplet_props);
      const double total_rate =
        inner_product(begin(triplet_props), end(triplet_props),
                      begin(the_model.triplet_rates), 0.0);
      cerr << "[expected mutations per site at root: "
           << total_rate << "]" << endl;
    }

    write_root_to_pathfile_global(pathfile, th.node_names[0], s);

    vector<state_seq> sequences(n_nodes, s);
    vector<size_t> events(8, 0);

    // iterate over the nodes in the tree
    for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
      const double curr_branch_len = th.branches[node_id];
      if (VERBOSE)
        cerr << "[SIMULATING: " << th.node_names[node_id]
             << " (" << curr_branch_len << ")]" << endl;

      TripletSampler ts(sequences[th.parent_ids[node_id]]);
      double time_value = 0;
      vector<GlobalJump> the_path;

      // sample changes along the current branch
      while (time_value < curr_branch_len)
        sample_jump(the_model, curr_branch_len, gen, ts, the_path, time_value,
                    events);

      // extract the sequence at the node
      ts.get_sequence(sequences[node_id]);

      append_to_pathfile_global(pathfile, th.node_names[node_id], the_path);

      if (VERBOSE)
        cerr << "[SUMMARY:]" << endl
             << summary_string(sequences[node_id]) << endl;
    }

    /* report frequencies of simulated events */
    if (VERBOSE) {
      cerr << "[FREQUENCIES OF SAMPLED EVENTS]" << endl;
      for (size_t i = 0; i < 8; ++i)
        cerr << std::bitset<3>(i) << '\t' << events[i] << endl;
      const double total_events = accumulate(begin(events), end(events), 0.0);
      cerr << "[TOTAL SAMPLED EVENTS: " << total_events << "]" << endl
           << "[EVENTS PER SITE PER UNIT TIME: " <<
        total_events/(evolutionary_time*n_sites) << "]" << endl;
    }

    if (VERBOSE)
      cerr << "[WRITING EPIGENOMIC STATES]" << endl;
    write_output(write_only_leaves, th, outfile, th.node_names, sequences);

  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
