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
#include "SingleSiteSampler.hpp"
#include "ContinuousTimeMarkovModel.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::min;
using std::runtime_error;

////////////////////////////////////////////////////////////////////////////////
//////////   copied/modified from WTF                                 //////////
////////////////////////////////////////////////////////////////////////////////
// collect rates and interval lengths
static void
collect_segment_info(const vector<double> &triplet_rates,
                     const Path &l, const Path &r,
                     vector<SegmentInfo> &seg_info) {
  Environment env(l, r);
  const size_t n_intervals = env.left.size();
  seg_info = vector<SegmentInfo>(n_intervals);

  for (size_t i = 0; i < n_intervals; ++i) {
    const size_t pattern0 = triple2idx(env.left[i], false, env.right[i]);
    const size_t pattern1 = triple2idx(env.left[i], true, env.right[i]);
    seg_info[i] = SegmentInfo(triplet_rates[pattern0], triplet_rates[pattern1],
                              env.breaks[i] - (i == 0 ? 0.0 : env.breaks[i-1]));
    assert(seg_info[i].len > 0.0);
  }
}

static double
root_post_prob0(const size_t site_id, const vector<Path> &the_paths,
                const vector<vector<double> > &horiz_tr_prob,
                const vector<double> &q) {

  const size_t left_st = the_paths[site_id - 1].init_state;
  const size_t right_st = the_paths[site_id + 1].init_state;

  const double p0 = (horiz_tr_prob[left_st][0]*horiz_tr_prob[0][right_st])*q[0];
  const double p1 = (horiz_tr_prob[left_st][1]*horiz_tr_prob[1][right_st])*q[1];

  return p0/(p0 + p1);
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
forward_sample_branch(const vector<SegmentInfo> &seg_info,
                      const bool par_state, bool &ch_state, std::mt19937 &gen,
                      vector<bool> &bp_states_branch) {

  const size_t n_intervals = seg_info.size();
  assert(bp_states_branch.size() == n_intervals);
  // propose a new path
  bool prev_state = par_state;
  ch_state = par_state;

  for (size_t m = 0; m < n_intervals; ++m) {
    const vector<double> rates {seg_info[m].rate0, seg_info[m].rate1};
    forward_sample_interval(rates, prev_state, seg_info[m].len, gen, ch_state);
    prev_state = ch_state;
    bp_states_branch[m] = ch_state;
  }
}


static void
downward_sampling_fs(const TreeHelper &th,
                     const vector<vector<SegmentInfo> > &seg_info,
                     const bool root_state, const vector<bool> &leaves_state,
                     std::mt19937 &gen, vector<vector<size_t> > &state_counts,
                     const size_t max_iterations) {

  const size_t n_nodes = th.subtree_sizes.size();
  bool success = false;
  size_t num_sampled = 0;

  vector<vector<bool> > bp_states;
  // initialize bp_states, entries will be rewritten during sampling
  for (size_t node_id = 0; node_id < n_nodes; ++node_id) {
    const size_t n_intervals = seg_info[node_id].size();
    vector<bool> bp_states_branch(n_intervals, false);
    bp_states.push_back(bp_states_branch);
  }

  while(!success && num_sampled < max_iterations) {
    vector<bool> nodes_state(n_nodes, false);
    nodes_state[0] = root_state;
    bool branch_end_state;

    bool leaves_sampled_ok = true;

    // preorder traversal of the tree
    for (size_t node_id = 1; node_id < n_nodes && leaves_sampled_ok;
         ++node_id) {
      forward_sample_branch(seg_info[node_id],
                            nodes_state[th.parent_ids[node_id]],
                            branch_end_state, gen, bp_states[node_id]);
      nodes_state[node_id] = branch_end_state;

      // check if leaf state sampled is valid
      if (is_leaf(th.subtree_sizes[node_id]))
          leaves_sampled_ok = (nodes_state[node_id] == leaves_state[node_id]);
    }
    ++num_sampled;
    success = leaves_sampled_ok;
  }

  // if success append jump_times to new_path
  if (success)
    for (size_t node_id = 0; node_id < n_nodes; ++node_id)
      for (size_t i = 0; i < seg_info[node_id].size(); ++i)
        state_counts[node_id][i] += !bp_states[node_id][i];
}


////////////////////////////////////////////////////////////////////////////////
//////////   Sample break points' state from posterior probabilities  //////////
////////////////////////////////////////////////////////////////////////////////
// posterior sampling
static void
posterior_sampling(const TreeHelper &th,
                   const vector<vector<SegmentInfo> > &seg_info,
                   const vector<FelsHelper> &fh, const bool root_state,
                   std::mt19937 &gen, vector<vector<size_t> > &state_counts) {

  const size_t n_nodes = th.subtree_sizes.size();

  for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
    const size_t n_intervals = seg_info[node_id].size();
    bool par_state = root_state;

    for (size_t m = 0; m < n_intervals; ++m) {
      // compute conditional posterior probability
      vector<vector<double> > P; // transition prob matrix
      continuous_time_trans_prob_mat(seg_info[node_id][m].rate0,
                                     seg_info[node_id][m].rate1,
                                     seg_info[node_id][m].len, P);

      const double p0 = P[par_state][0] / fh[node_id].p[m][par_state] *
      ((m == n_intervals - 1) ? fh[node_id].q[0] : fh[node_id].p[m + 1][0]);

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

    if (VERBOSE)
      cerr << "TEST SITE: " << site << endl
      << "PATHS TO SIMULATE: " << n_paths_to_sample << endl;

    if (VERBOSE)
      cerr << "[READING PATHS: " << pathsfile << "]" << endl;
    vector<vector<Path> > paths; // along multiple branches
    vector<string> node_names;
    read_paths(pathsfile, node_names, paths);
    assert(paths.back().size() == 5); // require 5 sites

    TreeHelper th;
    if (tree_file.empty()) {
      // test single edge
      if (VERBOSE)
        cerr << "[TREE NOT SPECIFIED: WILL LOAD FIRST PATH AS SINGLE BRANCH]"
             << endl;

      paths.resize(2);
      node_names.resize(2);

      th = TreeHelper(paths.back()[site].tot_time);
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
    assert(node_names.size() == th.n_nodes &&
           paths.size() == th.n_nodes);

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

    vector<vector<SegmentInfo> > seg_info(n_nodes);
    for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
      collect_segment_info(the_model.triplet_rates,
                           paths[node_id][site - 1],
                           paths[node_id][site + 1], seg_info[node_id]);
    }

    // (1) Upward pruning
    vector<FelsHelper> fh;
    pruning(th, site, paths, seg_info, fh);

    // counts of mid state 0 at breakpoints n_branches * n_intervals
    vector<vector<size_t> > all_bp_state0_fs, all_bp_state0;
    vector<vector<size_t> > all_bp_state1_fs, all_bp_state1;
    for (size_t node_id = 0; node_id < n_nodes; ++node_id) {
      const size_t n_intervals = seg_info[node_id].size();
      vector<size_t> bp_state0_fs(n_intervals, 0), bp_state0(n_intervals, 0);
      vector<size_t> bp_state1_fs(n_intervals, 0), bp_state1(n_intervals, 0);
      all_bp_state0_fs.push_back(bp_state0_fs);
      all_bp_state0.push_back(bp_state0);
      all_bp_state1_fs.push_back(bp_state1_fs);
      all_bp_state1.push_back(bp_state1);
    }

    // collect leaf nodes' state
    vector<bool> leaves_state(th.subtree_sizes.size(), false);
    for (size_t node_id = 0; node_id < th.subtree_sizes.size(); ++node_id)
      if (is_leaf(th.subtree_sizes[node_id]))
        leaves_state[node_id] = paths[node_id][site].end_state();

    size_t n_paths_from_zero = 0, n_paths_sampled = 0;

    while (n_paths_sampled < n_paths_to_sample) {
      if (VERBOSE && n_paths_sampled * 10 % n_paths_to_sample == 0)
        cerr << "FINISHED: " << n_paths_sampled * 100 / n_paths_to_sample
             << '%' << endl;

      // sample new root state

      const double root_p0 = root_post_prob0(site, paths[1], the_model.init_T,
                                             fh[0].q);
      std::uniform_real_distribution<double> unif(0.0, 1.0);
      bool new_root_state = (unif(gen) > root_p0);

      // Recursively sample the whole tree
      if (new_root_state) { // root state is 1
        downward_sampling_fs(th, seg_info, new_root_state, leaves_state,
                             gen, all_bp_state1_fs, max_iterations);

        posterior_sampling(th, seg_info, fh, new_root_state, gen,
                           all_bp_state1);
      } else { // root state is 0
        downward_sampling_fs(th, seg_info, new_root_state, leaves_state,
                             gen, all_bp_state0_fs, max_iterations);

        posterior_sampling(th, seg_info, fh, new_root_state, gen,
                           all_bp_state0);
        ++n_paths_from_zero;
      }

      ++n_paths_sampled;
    }
    cerr << "FINISHED: " << n_paths_sampled * 100 / n_paths_to_sample
         << '%' << endl;

    // write output
    std::ofstream out(outfile.c_str());
    out << "NODE" << '\t' << "BREAK" << '\t'
    << "LEFT_INIT_STATE" << '\t' << "RIGHT_INIT_STATE" << '\t'
    << "START0_MID_0_PROB_FS" << '\t'
    << "START0_MID_0_PROB_PR" << '\t'
    << "START1_MID_1_PROB_FS" << '\t'
    << "START1_MID_1_PROB_PR" << endl;

    for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
      Environment env(paths[node_id][site-1], paths[node_id][site+1]);
      for (size_t i = 0; i < seg_info[node_id].size(); ++i) {
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
