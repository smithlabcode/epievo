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

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

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

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::runtime_error;
using std::numeric_limits;

static void
assign_changes_to_sites(const vector<GlobalJump> &global_path,
                        vector<Path> &by_site) {

  const size_t n_changes = global_path.size();
  for (size_t i = 0; i < n_changes; ++i)
    by_site[global_path[i].position].jumps.push_back(global_path[i].timepoint);
}


///////////////////////////////////////////////////////////////////////////

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


int main(int argc, const char **argv) {
  try {

    static const size_t max_sample_count = 10000;

    bool VERBOSE = false;
    string outfile;

    size_t n_paths_to_sample = 1000;

    size_t n_sites = 5;
    double evo_time = 1.0;

    string target("00000");

    size_t rng_seed = std::numeric_limits<size_t>::max();

    ///////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]), "test end-conditioned samplers",
                           "<param-file>");
    opt_parse.add_opt("sites", 's', "number of sites", false, n_sites);
    opt_parse.add_opt("time", 't', "duration of time interval", false, evo_time);
    opt_parse.add_opt("leaf", 'l', "target leaf sequence", false,
                      target);
    opt_parse.add_opt("samples", 'n', "number of paths to sample",
                      false, n_paths_to_sample);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("seed", 'S', "rng seed", false, rng_seed);
    opt_parse.add_opt("outfile", 'o', "outfile (sampling statistics)",
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    const string param_file(leftover_args.front());
    ///////////////////////////////////////////////////////////////////////////

    /* (1) INITIALIZING THE (FAKE) TREE */
    if (VERBOSE)
      cerr << "[INITIALIZING TWO NODE TREE (time= "
           << evo_time << ")]" << endl;
    const TreeHelper th(evo_time);

    /* (2) READING PARAMETERS FROM FILE */
    if (VERBOSE)
      cerr << "[READING PARAMETER FILE: " << param_file << "]" << endl;
    EpiEvoModel the_model;
    read_model(param_file, the_model);
    // the_model.scale_triplet_rates();
    if (VERBOSE)
      cerr << the_model << endl;

    /* (3) INITIALIZE THE RANDOM NUMBER GENERATOR */
    /* standard mersenne_twister_engine seeded with rd() */
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    std::mt19937 gen(rng_seed);
    if (VERBOSE)
      cerr << "[INITIALIZED RNG (seed=" << rng_seed << ")]" << endl;

    if (VERBOSE)
      cerr << "[SAMPLING HISTORY FOR FULL SEQUENCE]" << endl;
    if (VERBOSE)
      cerr << "[SIMULATING ROOT: " << th.node_names[0]
           << " (N SITES=" << n_sites << ")]" << endl;
    vector<char> root_seq = {false, false, false, false, false};
    n_sites = 5;
    size_t the_site = 2;
    // the_model.sample_state_sequence_init(n_sites, gen, root_seq);
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

    StateSeq target_leaf_seq(s);

    cerr << "target=" << target_leaf_seq << endl;

    /* (4) SAMPLE THE PATHS GLOBALLY */
    size_t target_count = 0;
    for (size_t i = 0; i < n_paths_to_sample; ++i) {
      TripletSampler ts(s);

      vector<GlobalJump> the_path;
      double time_value = 0;

      while (time_value < evo_time)
        sample_jump(the_model, evo_time, gen, ts, the_path, time_value);

      StateSeq leaf_seq;
      ts.get_sequence(leaf_seq);

      // check if we have the desired leaf state
      if (leaf_seq.seq == target_leaf_seq.seq) {
        ++target_count;

        vector<Path> path_by_site(n_sites);
        for (size_t j = 0; j < path_by_site.size(); ++j) {
          path_by_site[j].init_state = s.seq[j];
          path_by_site[j].tot_time = th.branches[1];
        }
        assign_changes_to_sites(the_path, path_by_site);

        vector<vector<Path> > paths(2, path_by_site);
        vector<Path> proposed_path;
        Gibbs_site(the_model, th, the_site, paths, gen, proposed_path);
      }
    }
    if (VERBOSE)
      cerr << "success rate in forward simulation: "
           << static_cast<double>(target_count)/n_paths_to_sample << endl;
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
