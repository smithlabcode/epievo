/* Copyright (C) 2018 University of Southern California
 *                    Jianghan Qu and Andrew D Smith
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

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "PhyloTreePreorder.hpp"
#include "Path.hpp"  /* related to Path */
#include "EpiEvoModel.hpp" /* model_param */
#include "StateSeq.hpp"
#include "SingleSampler.hpp"
#include "ContinuousTimeMarkovModel.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::min;
using std::runtime_error;

////////////////////////////////////////////////////////////////////////////////
//////////   copied/modified from SingleSampler.cpp                   //////////
////////////////////////////////////////////////////////////////////////////////
// collect rates and interval lengths
static void
rates_on_branch(const vector<double> &triplet_rates,
                const Path &left_neighbor_path,
                const Path &right_neighbor_path,
                vector<vector<double> > &interval_rates,
                vector<double> &interval_lengths) {

  const Environment env(left_neighbor_path, right_neighbor_path);

  const size_t n_intervals = env.left.size();

  interval_rates = vector<vector<double> > (n_intervals, vector<double>(2, 0.0));
  interval_lengths = vector<double>(n_intervals, 0.0);

  for (size_t i = 0; i < n_intervals; ++i) {

    const size_t pattern0 = triple2idx(env.left[i], 0ul, env.right[i]);
    const size_t pattern1 = flip_mid_bit(pattern0);

    interval_rates[i][0] = triplet_rates[pattern0];
    interval_rates[i][1] = triplet_rates[pattern1];

    interval_lengths[i] = (i == 0) ?
      env.breaks[0] : env.breaks[i] - env.breaks[i - 1];

    assert(interval_lengths[i] > 0);
  }
}


static double
root_post_prob0(const size_t site,
                const vector<vector<Path> > &all_paths,
                const vector<vector<double> > &root_trans_prob,
                const vector<vector<vector<double> > > &all_p) {
  // compute posterior probability at root node
  // (i.e. the init_state of all children)
  size_t lstate = all_paths[1][site - 1].init_state;
  size_t rstate = all_paths[1][site + 1].init_state;
  double p0 = root_trans_prob[lstate][0] * root_trans_prob[0][rstate];
  double p1 = root_trans_prob[lstate][1] * root_trans_prob[1][rstate];
  p0 *= all_p[1][0][0];
  p1 *= all_p[1][0][1];

  double root_p0 = p0 / (p0+p1);
  return root_p0;
}


////////////////////////////////////////////////////////////////////////////////
//////////   Downward_sampling_branch Forward sampling                //////////
////////////////////////////////////////////////////////////////////////////////
static void
forward_sample_interval(const vector<double> &rates,
                        const bool initial_state, const double tot_time,
                        std::mt19937 &gen, bool &end_state) {

  double time_value = 0;

  end_state = initial_state;

  while (time_value < tot_time) {
    std::exponential_distribution<double> exp_distr(rates[end_state]);
    const double holding_time =
      std::max(exp_distr(gen), std::numeric_limits<double>::min());

    time_value += holding_time;

    if (time_value < tot_time) {
      end_state = complement_state(end_state);
    }
  }
}


// using Forward sampling + rejection
static void
downward_sampling_branch_fs(const vector<vector<double> > &interval_rates,
                            const vector<double> &interval_lengths,
                            const bool is, const bool leaf_state,
                            std::mt19937 &gen,
                            vector<size_t> &state_counts,
                            const size_t max_iterations) {

  const size_t n_intervals = interval_rates.size();
  bool reach_leaf_state = false;
  size_t num_sampled = 0;
  vector<bool> state_proposed;

  while(!reach_leaf_state && num_sampled < max_iterations) {
    state_proposed.clear();
    // propose a new path
    bool par_state = is;
    bool new_state;

    double time_passed = 0;

    for (size_t m = 0; m < n_intervals; ++m) {
      forward_sample_interval(interval_rates[m], par_state,
                              interval_lengths[m], gen, new_state);
      time_passed += interval_lengths[m];
      par_state = new_state;
      state_proposed.push_back(new_state);
    }
    ++num_sampled;
    reach_leaf_state = (new_state == leaf_state);
  }

  // if success append jump_times to new_path
  if (reach_leaf_state)
    for (size_t i = 0; i < n_intervals; ++i)
      state_counts[i] += !state_proposed[i];
}



// posterior sampling
static void
posterior_sampling(const vector<vector<double> > &interval_rates,
                   const vector<double> &interval_lengths,
                   const vector<vector<double> > &all_p,
                   const bool is, const bool leaf_state,
                   std::mt19937 &gen,
                   vector<size_t> &state_counts) {

  const size_t n_intervals = interval_rates.size();
  size_t par_state = is;

  for (size_t m = 0; m < n_intervals; ++m) {
    // compute conditional posterior probability
    vector<vector<double> > P; // transition prob matrix
    const CTMarkovModel ctmm(interval_rates[m]);
    ctmm.get_trans_prob_mat(interval_lengths[m], P);
    double p0_post = m < n_intervals - 1 ? all_p[m+1][0] : leaf_state == false;
    double p0 = p0_post*P[par_state][0]/all_p[m][par_state];

    // generate random state at break point
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    bool new_state = (unif(gen) > p0);
    if (!new_state)
      state_counts[m]++;

    par_state = new_state;
  }
}


////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char **argv) {
  try {

    bool VERBOSE = false;
    string outfile;
    string outstatefile;

    size_t max_iterations = 1000000;
    size_t n_paths_to_sample = 1000;
    size_t test_site = 2; // we test 5-site path
    size_t test_branch = 1; // branch to test in the tree

    string param_file;
    string tree_file;

    size_t rng_seed = std::numeric_limits<size_t>::max();

    ///////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]), "test triple path",
                           "<paths-file>");
    opt_parse.add_opt("param", 'p', "parameter file",
                      true, param_file);
    opt_parse.add_opt("tree", 't', "tree file in newick format",
                      true, tree_file);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);
    opt_parse.add_opt("paths", 'n', "number of paths to sample",
                      false, n_paths_to_sample);
    opt_parse.add_opt("outfile", 'o', "outfile (sampling statistics)",
                      true, outfile);
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
    const string pathsfile(leftover_args.front());
    ///////////////////////////////////////////////////////////////////////////

    if (VERBOSE)
      cerr << "[READING PARAMETER FILE: " << param_file << ", "
           << tree_file << "]" << endl;

    EpiEvoModel the_model;
    read_model(param_file, the_model);
    if (VERBOSE)
      cerr << the_model << endl;
    
    if (VERBOSE)
      cerr << "[READING TREE FILE: " << tree_file << "]" << endl;
    PhyloTreePreorder the_tree; // tree topology and branch lengths
    std::ifstream tree_in(tree_file.c_str());
    if (!tree_in || !(tree_in >> the_tree))
      throw std::runtime_error("bad tree file: " + tree_file);
    const size_t n_nodes = the_tree.get_size();
    TreeHelper th(the_tree);
    // remove undesired branch
    th.subtree_sizes.erase(th.subtree_sizes.begin()+2);
    th.node_names.erase(th.node_names.begin()+2);
    th.parent_ids.erase(th.parent_ids.begin()+2);
    th.branches.erase(th.branches.begin()+2);


    if (VERBOSE)
      cerr << "[READING PATHS: " << pathsfile << "]" << endl;
    vector<vector<Path> > all_paths; // along multiple branches
    vector<string> node_names;
    read_paths(pathsfile, node_names, all_paths);
    // remove undesired branch
    all_paths.erase(all_paths.begin()+2);
    node_names.erase(node_names.begin()+2);

    if (VERBOSE)
      cerr << "TEST BRANCH: " << test_branch << endl
           << "TEST SITE: " << test_site << endl
           << "PATHS TO SIMULATE: " << n_paths_to_sample << endl;
    // standard mersenne_twister_engine seeded with rd()
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    if (VERBOSE)
      cerr << "rng seed: " << rng_seed << endl;

    std::mt19937 gen(rng_seed);

    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    cerr << "----- TEST upward downward sampling BELOW ---------" << endl;

    vector<vector<vector<double> > > all_interval_rates(n_nodes);
    vector<vector<double> > all_interval_lengths(n_nodes);
    rates_on_branch(the_model.triplet_rates,
                    all_paths[1][test_site - 1],
                    all_paths[1][test_site + 1],
                    all_interval_rates[1],
                    all_interval_lengths[1]);

    // (1) Upward pruning
    vector<vector<vector<double> > > all_p;
    all_p.resize(th.subtree_sizes.size());
    pruning(the_model.triplet_rates, th.subtree_sizes, test_site,
            all_paths, all_interval_rates, all_interval_lengths, all_p);

    // counts of mid state 0 at breakpoints
    vector<size_t> bp_state0_fs, bp_state0;
    bp_state0_fs.resize(all_interval_rates[1].size(), 0);
    bp_state0.resize(all_interval_rates[1].size(), 0);

    size_t n_paths_sampled = 0;

    while (n_paths_sampled < n_paths_to_sample) {
      if (VERBOSE && n_paths_sampled * 10 % n_paths_to_sample == 0)
        cerr << "FINISHED: " << n_paths_sampled * 100 / n_paths_to_sample
             << '%' << endl;

      // (2) Foward sample a path conditioned on given leaf state
      Path new_path_fs;

      // sample new root state
      const double root_p0 = root_post_prob0(test_site, all_paths,
                                             the_model.init_T, all_p);
      std::uniform_real_distribution<double> unif(0.0, 1.0);
      bool new_root_state = (unif(gen) > root_p0);

      new_path_fs.init_state = new_root_state;
      new_path_fs.tot_time = all_paths[1][test_site - 1].tot_time;

      downward_sampling_branch_fs(all_interval_rates[1],
                                  all_interval_lengths[1],
                                  new_root_state,
                                  all_paths[1][test_site].end_state(),
                                  gen, bp_state0_fs, max_iterations);

      // (3) Posterior sampling:
      posterior_sampling(all_interval_rates[1],
                         all_interval_lengths[1], all_p[1],
                         new_root_state,
                         all_paths[1][test_site].end_state(),
                         gen, bp_state0);

      ++n_paths_sampled;
    }
    cerr << "FINISHED: " << n_paths_sampled * 100 / n_paths_to_sample
                         << '%' << endl;

    // write output
    std::ofstream out(outfile.c_str());
    out << "BREAK" << '\t'
        << "LEFT_INIT_STATE" << '\t' << "RIGHT_INIT_STATE" << '\t'
        << "PROP_MID_END_0_FS" << '\t' << "PROP_MID_END_0_PR" << endl;

    Environment env(all_paths[1][test_site-1], all_paths[1][test_site+1]);
    for (size_t i = 0; i < all_interval_lengths[1].size(); ++i) {
      out << i+1 << '\t'
          << env.left[i] << '\t' << env.right[i] << '\t'
          << 1.0 * bp_state0_fs[i] / n_paths_sampled << '\t'
          << 1.0 * bp_state0[i] / n_paths_sampled << endl;
    }
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
