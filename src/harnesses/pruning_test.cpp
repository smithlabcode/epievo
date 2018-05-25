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
#include "TreeHelper.hpp"
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
root_post_prob0(const vector<size_t> &children,
                const size_t site,
                const vector<vector<Path> > &all_paths,
                const vector<vector<double> > &root_trans_prob,
                const vector<vector<vector<double> > > &all_p) {
  // compute posterior probability at root node
  // (i.e. the init_state of all children)
  size_t lstate = all_paths[children[0]][site - 1].init_state;
  size_t rstate = all_paths[children[0]][site + 1].init_state;
  double p0 = root_trans_prob[lstate][0] * root_trans_prob[0][rstate];
  double p1 = root_trans_prob[lstate][1] * root_trans_prob[1][rstate];
  for (size_t idx = 0; idx < children.size(); ++idx) {
    p0 *= all_p[children[idx]][0][0];
    p1 *= all_p[children[idx]][0][1];
  }
  
  double root_p0 = p0 / (p0+p1);
  return root_p0;
}


////////////////////////////////////////////////////////////////////////////////
//////////   Forward sample paths conditioned on leaves' state        //////////
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


// BFS forward sampling
// return fail if propose wrong leaf state
static void
forward_sample_branch(const vector<vector<double>> &interval_rates,
                      const vector<double> &interval_lengths,
                      const bool par_state, bool &ch_state, std::mt19937 &gen,
                      vector<bool> &bp_states_branch) {

  assert(bp_states_branch.size() == interval_lengths.size());
  // propose a new path
  bool prev_state = par_state;
  ch_state = par_state;
  
  const size_t n_intervals = interval_lengths.size();
  
  for (size_t m = 0; m < n_intervals; ++m) {
    forward_sample_interval(interval_rates[m], prev_state,
                            interval_lengths[m], gen, ch_state);
    prev_state = ch_state;
    bp_states_branch[m] = ch_state;
  }
}


static void
downward_sampling_fs(const vector<vector<vector<double>>> &all_interval_rates,
                     const vector<vector<double>> &all_interval_lengths,
                     const vector<size_t> &subtree_sizes,
                     const vector<size_t> &parent_ids,
                     const bool root_state, const vector<bool> &leaves_state,
                     std::mt19937 &gen, vector<vector<size_t> > &state_counts,
                     const size_t max_iterations) {
  
  const size_t n_nodes = subtree_sizes.size();
  bool success = false;
  size_t num_sampled = 0;
  
  vector<vector<bool> > bp_states;
  // initialize bp_states, entries will be rewritten during sampling
  for (size_t node_id = 0; node_id < n_nodes; ++node_id) {
    const size_t n_intervals = all_interval_lengths[node_id].size();
    vector<bool> bp_states_branch(n_intervals, false);
    bp_states.push_back(bp_states_branch);
  }
  
  while(!success && num_sampled < max_iterations) {
    vector<bool> nodes_state(n_nodes, false);
    nodes_state[0] = root_state;
    
    bool leaves_sampled_ok = false;
    
    // preorder traversal of the tree
    for (size_t node_id = 1; node_id < n_nodes && leaves_sampled_ok;
         ++node_id) {
      forward_sample_branch(all_interval_rates[node_id],
                            all_interval_lengths[node_id],
                            nodes_state[parent_ids[node_id]],
                            nodes_state[node_id], gen, bp_states[node_id]);
      
      // check if leaf state sampled is valid
      if (is_leaf(subtree_sizes[node_id]))
          leaves_sampled_ok = (nodes_state[node_id] == leaves_state[node_id]);
    }
    ++num_sampled;
    success = leaves_sampled_ok;
  }
  
  // if success append jump_times to new_path
  if (success)
    for (size_t node_id = 0; node_id < n_nodes; ++node_id)
      for (size_t i = 0; i < all_interval_lengths[node_id].size(); ++i)
        state_counts[node_id][i] += !bp_states[node_id][i];
}


////////////////////////////////////////////////////////////////////////////////
//////////   Sample break points' state from posterior probabilities  //////////
////////////////////////////////////////////////////////////////////////////////
// posterior sampling
static void
posterior_sampling(const vector<vector<vector<double>>> &all_interval_rates,
                   const vector<vector<double>> &all_interval_lengths,
                   const vector<vector<vector<double> > > &all_p,
                   const vector<size_t> &subtree_sizes,
                   const vector<size_t> &parent_ids,
                   const bool root_state, const vector<bool> &leaves_state,
                   std::mt19937 &gen, vector<vector<size_t> > &state_counts) {

  const size_t n_nodes = subtree_sizes.size();

  for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
    const size_t n_intervals = all_interval_lengths[node_id].size();
    bool par_state = root_state;
    
    for (size_t m = 0; m < n_intervals; ++m) {
      // compute conditional posterior probability
      vector<vector<double> > P; // transition prob matrix
      const CTMarkovModel ctmm(all_interval_rates[node_id][m]);
      ctmm.get_trans_prob_mat(all_interval_lengths[node_id][m], P);
      
      double p0_post = 1.0;
      if (m == (n_intervals - 1))
        if (is_leaf(subtree_sizes[node_id])) {
          vector<size_t> children;
          get_children(node_id, subtree_sizes, children);
          for (size_t idx = 0; idx < children.size(); ++idx)
            p0_post *= all_p[children[idx]][0][0];
        }
        else
          p0_post = (leaves_state[node_id] == false);
      else
        p0_post = all_p[node_id][m+1][0];
      
      double p0 = p0_post*P[par_state][0]/all_p[node_id][m][par_state];
      
      // generate random state at break point
      std::uniform_real_distribution<double> unif(0.0, 1.0);
      bool new_state = (unif(gen) > p0);
      if (!new_state)
        state_counts[node_id][m]++;
      
      par_state = new_state;
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char **argv) {
  try {

    static const size_t site = 2;

    bool VERBOSE = false;
    string outfile;
    string outstatefile;

    size_t max_iterations = 1000000;
    size_t n_paths_to_sample = 1000;
    size_t site = 2; // we test 5-site path

    string param_file;
    string tree_file;

    size_t rng_seed = std::numeric_limits<size_t>::max();

    ///////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]), "test triple path",
                           "<paths-file>");
    opt_parse.add_opt("param", 'p', "parameter file",
                      true, param_file);
    opt_parse.add_opt("tree", 't', "tree file in newick format",
                      false, tree_file);
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
      cerr << "[READING PARAMETER FILE: " << param_file << "]" << endl;

    EpiEvoModel the_model;
    read_model(param_file, the_model);
    if (VERBOSE)
      cerr << the_model << endl;
    
    vector<vector<Path> > all_paths; // n_nodes * n_intervals
    vector<string> node_names;
    read_paths(pathsfile, node_names, all_paths);
    
    TreeHelper th;
    
    if (tree_file.empty()) {
      // test single edge
      if (VERBOSE)
        cerr << "[TREE NOT SPECIFIED: WILL LOAD FIRST PATH AS SINGLE BRANCH]"
             << endl;
      
      all_paths.resize(2);
      node_names.resize(2);
      
      th = TreeHelper(all_paths.back()[site].tot_time);
      th.node_names = node_names;
    } else {
      // test whole tree
      if (VERBOSE)
        cerr << "[READING TREE FILE: " << tree_file << "]" << endl;
      std::ifstream tree_in(tree_file.c_str());
      PhyloTreePreorder the_tree;
      if (!tree_in || !(tree_in >> the_tree))
        throw std::runtime_error("bad tree file: " + tree_file);
      th = TreeHelper(the_tree);
    }

    if (VERBOSE)
      cerr << "TEST SITE: " << site << endl
           << "PATHS TO SIMULATE: " << n_paths_to_sample << endl;

    if (VERBOSE)
      cerr << "[READING TREE FILE: " << tree_file << "]" << endl;
    PhyloTreePreorder the_tree; // tree topology and branch lengths
    std::ifstream tree_in(tree_file.c_str());
    if (!tree_in || !(tree_in >> the_tree))
      throw std::runtime_error("bad tree file: " + tree_file);
    TreeHelper th(the_tree);
    assert(th.n_nodes == 3ul);

    if (VERBOSE)
      cerr << "[READING PATHS: " << pathsfile << "]" << endl;
    vector<vector<Path> > all_paths; // along multiple branches
    vector<string> node_names;
    read_paths(pathsfile, node_names, all_paths);
    assert(node_names.size() == th.n_nodes &&
           all_paths.size() == th.n_nodes);
    assert(all_paths.back().size() == 5); // require 5 sites

    vector<Path> the_paths(all_paths[1]);

    if (VERBOSE)
      cerr << "number of paths to simulate: " << n_paths_to_sample << endl;

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
    
    const size_t n_nodes = th.n_nodes;

    vector<vector<vector<double> > > all_interval_rates(th.n_nodes);
    vector<vector<double> > all_interval_lengths(th.n_nodes);
    for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
      rates_on_branch(the_model.triplet_rates,
                      all_paths[node_id][site - 1],
                      all_paths[node_id][site + 1],
                      all_interval_rates[node_id],
                      all_interval_lengths[node_id]);
    }

    // (1) Upward pruning
    vector<vector<vector<double> > > all_p;
    all_p.resize(th.subtree_sizes.size());
    pruning(the_model.triplet_rates, th.subtree_sizes, site,
            all_paths, all_interval_rates, all_interval_lengths, all_p);
    
    // counts of mid state 0 at breakpoints n_branches * n_intervals
    vector<vector<size_t> > all_bp_state0_fs, all_bp_state0;
    vector<vector<size_t> > all_bp_state1_fs, all_bp_state1;
    for (size_t node_id = 0; node_id < n_nodes; ++node_id) {
      const size_t n_intervals = all_interval_lengths[node_id].size();
      vector<size_t> bp_state0_fs(n_intervals, 0), bp_state0(n_intervals, 0);
      vector<size_t> bp_state1_fs(n_intervals, 0), bp_state1(n_intervals, 0);
      all_bp_state0_fs.push_back(bp_state0_fs);
      all_bp_state0.push_back(bp_state0);
      all_bp_state1_fs.push_back(bp_state1_fs);
      all_bp_state1.push_back(bp_state1);
    }
    
    size_t n_paths_from_zero = 0, n_paths_sampled = 0;
    
    while (n_paths_sampled < n_paths_to_sample) {
      if (VERBOSE && n_paths_sampled * 10 % n_paths_to_sample == 0)
        cerr << "FINISHED: " << n_paths_sampled * 100 / n_paths_to_sample
        << '%' << endl;
      
      // sample new root state
      const size_t root_id = 0; //all_paths[0] is empty
      vector<size_t> children;
      get_children(root_id, th.subtree_sizes, children);
      
      const double root_p0 = root_post_prob0(children, site, all_paths,
                                             the_model.init_T, all_p);
      std::uniform_real_distribution<double> unif(0.0, 1.0);
      bool new_root_state = (unif(gen) > root_p0);
      
      // collect leaf nodes' state
      vector<bool> leaves_state(th.subtree_sizes.size(), false);
      for (size_t node_id = 0; node_id < th.subtree_sizes.size(); ++node_id)
        if (is_leaf(th.subtree_sizes[node_id]))
          leaves_state[node_id] = all_paths[node_id][site].end_state();
      
      // Recursively sample the whole tree
      if (new_root_state) { // root state is 1
        downward_sampling_fs(all_interval_rates, all_interval_lengths,
                             th.subtree_sizes, th.parent_ids,
                             new_root_state, leaves_state,
                             gen, all_bp_state1_fs, max_iterations);
        
        posterior_sampling(all_interval_rates, all_interval_lengths, all_p,
                           th.subtree_sizes, th.parent_ids,
                           new_root_state, leaves_state,
                           gen, all_bp_state1);
      } else { // root state is 0
        downward_sampling_fs(all_interval_rates, all_interval_lengths,
                             th.subtree_sizes, th.parent_ids,
                             new_root_state, leaves_state,
                             gen, all_bp_state0_fs, max_iterations);
        
        posterior_sampling(all_interval_rates, all_interval_lengths, all_p,
                           th.subtree_sizes, th.parent_ids,
                           new_root_state, leaves_state,
                           gen, all_bp_state0);
        ++n_paths_from_zero;
      }
      
      ++n_paths_sampled;
      
      cerr << "FINISHED: " << n_paths_sampled * 100 / n_paths_to_sample
      << '%' << endl;
    }
    
    // write output
    std::ofstream out(outfile.c_str());
    for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
      out << "NODE" << '\t' << "BREAK" << '\t'
          << "LEFT_INIT_STATE" << '\t' << "RIGHT_INIT_STATE" << '\t'
          << "START0_MID_0_PROB_FS" << '\t'
          << "START0_MID_0_PROB_PR" << '\t'
          << "START1_MID_1_PROB_FS" << '\t'
          << "START1_MID_1_PROB_PR" << endl;
      
      Environment env(all_paths[node_id][site-1], all_paths[node_id][site+1]);
      for (size_t i = 0; i < all_interval_lengths[node_id].size(); ++i) {
        out << th.node_names[node_id] << '\t' << i+1 << '\t'
            << env.left[i] << '\t' << env.right[i] << '\t' << "0\t"
            << 1.0 * all_bp_state0_fs[node_id][i] / n_paths_from_zero << '\t'
            << 1.0 * all_bp_state0[node_id][i] / n_paths_from_zero << '\t'
            << 1.0 * all_bp_state1_fs[node_id][i] /
                       (n_paths_sampled - n_paths_from_zero)
            << '\t' << 1.0 * all_bp_state1[node_id][i] /
                               (n_paths_sampled - n_paths_from_zero)
            << endl;
      }
    }
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
