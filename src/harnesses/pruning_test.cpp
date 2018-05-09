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
#include <bitset>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_randist.h> /* chi-squared test */
#include <gsl/gsl_cdf.h>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "PhyloTreePreorder.hpp"
#include "Path.hpp"  /* related to Path */
#include "EpiEvoModel.hpp" /* model_param */
#include "StateSeq.hpp"
#include "SingleSampler.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::min;
using std::runtime_error;
using std::bitset;

static const size_t N_TRIPLETS = 8;
static const double MINWAIT = 1e-8;

////////////////////////////////////////////////////////////////////////////////
//////////   Record summary statistics                                //////////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//////////   copied/modified from SingleSampler.cpp                   //////////
////////////////////////////////////////////////////////////////////////////////
// collect rates and interval lengths
static void
rates_on_branch(const vector<double> &triplet_rates,
                const Path &l, const Path &r,
                vector<vector<double> > &interval_rates,
                vector<double> &interval_lengths) {
  Environment env(l, r);
  const size_t n_intervals = env.left.size();
  interval_rates = vector<vector<double> > (n_intervals, vector<double>(2, 0.0));
  interval_lengths = vector<double>(n_intervals, 0.0);
  
  for (size_t i = 0; i < n_intervals; ++i) {
    const size_t pattern0 = 4 * (size_t)(env.left[i]) + (size_t)(env.right[i]);
    const size_t pattern1 = pattern0 + 2;
    interval_rates[i][0] = triplet_rates[pattern0];
    interval_rates[i][1] = triplet_rates[pattern1];
    interval_lengths[i] = (i == 0) ? env.breaks[0] : env.breaks[i] - env.breaks[i-1];
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
                        const size_t is, size_t &es,
                        const double tot_time, const double time_passed,
                        std::mt19937 &gen,
                        vector<double> &jump_times) {
  
  double time_value = 0;
  size_t curr_state = is;
  
  while (time_value < tot_time) {
    const double holding_rate = rates[curr_state];
    std::exponential_distribution<double> exp_distr(holding_rate);
    const double holding_time =
    std::max(exp_distr(gen), std::numeric_limits<double>::min());
    
    time_value += holding_time;
    
    if (time_value < tot_time) {
      jump_times.push_back(time_passed + time_value);
      curr_state = complement_state(curr_state);
    }
  }
  
  es = curr_state;
}


// using Forward sampling + rejection
static void
downward_sampling_branch_fs(const vector<vector<double> > &interval_rates,
                            const vector<double> &interval_lengths,
                            const size_t site,
                            const vector<Path> &paths,
                            const vector<vector<double> > &all_p,
                            std::mt19937 &gen,
                            Path &new_path,
                            const size_t max_iterations) {
  
  const size_t n_intervals = interval_rates.size();
  const size_t leaf_state = paths[site].end_state();
  bool reach_leaf_state = false;
  size_t num_sampled = 0;
  vector<double> jump_times;

  while(!reach_leaf_state && num_sampled < max_iterations) {
    // propose a new path
    jump_times.clear();
    size_t par_state = new_path.init_state;
    size_t new_state;

    double time_passed = 0;
    
    for (size_t m = 0; m < n_intervals; ++m) {
      forward_sample_interval(interval_rates[m],
                              par_state, new_state,
                              interval_lengths[m], time_passed,
                              gen, jump_times);
      time_passed += interval_lengths[m];
      par_state = new_state;
    }
    
    ++num_sampled;
    reach_leaf_state = new_state == leaf_state;
  }
  
  // if success append jump_times to new_path
  if (reach_leaf_state)
    for (size_t i = 0; i < jump_times.size(); ++i) {
      new_path.jumps.push_back(jump_times[i]) ;
    }
}


/*
// using Forward sampling + rejection
static void
downward_sampling_branch_fs(const vector<vector<double> > &interval_rates,
                            const vector<double> &interval_lengths,
                            const size_t site,
                            const vector<Path> &paths,
                            const vector<vector<double> > &all_p,
                            std::mt19937 &gen,
                            Path &new_path,
                            const size_t max_iterations) {
  
  const size_t n_intervals = interval_rates.size();
  const size_t leaf_state = paths[site].end_state();
  
  bool reach_target = false;
  
  // propose a new path
  vector<double> jump_times;
  size_t num_sampled = 0;
  
  size_t par_state = new_path.init_state;
  
  double time_passed = 0;
  
  for (size_t m = 0; m < n_intervals - 1; ++m) {
    // compute conditional posterior probability
    vector<vector<double> > P; // transition prob matrix
    trans_prob_mat(interval_rates[m][0], interval_rates[m][1],
                   interval_lengths[m], P);
    double p0 = (all_p[m+1][0] * P[par_state][0] /
                 all_p[m][par_state]);
    
    // generate random state at break point
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    bool new_state = (unif(gen) > p0);
    
    // generate path
    vector<double> jump_times;
    forward_sample(interval_rates[m],
                   par_state, (size_t)(new_state),
                   interval_lengths[m], gen, jump_times, max_iterations);
    
    // append jump_times to new_path
    for (size_t i = 0; i < jump_times.size(); ++i) {
      new_path.jumps.push_back(time_passed + jump_times[i]) ;
    }
    // prepare for next interval
    time_passed += interval_lengths[m];
    par_state = new_state;
  }
  
  // generate path
  vector<double> jump_times;
  const size_t m = n_intervals - 1;
  // prepare helper values
  forward_sample(interval_rates[m], par_state, leaf_state,
                 interval_lengths[m], gen, jump_times, max_iterations);
  // append jump_times to new_path
  for (size_t i = 0; i < jump_times.size(); ++i) {
    new_path.jumps.push_back(time_passed + jump_times[i]) ;
  }
  
  assert(new_path.end_state() == leaf_state);
  for (size_t i = 0; i < fs_jump_times.size(); ++i) {
    jump_times.push_back(fs_jump_times[i]);
  }
  
}
*/
////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char **argv) {
  try {

    bool VERBOSE = false;
    bool SCALE = false;
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
                           " <paths-file>");
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
    read_model(SCALE, param_file, tree_file, the_model);
    // remove undesired branch
    the_model.subtree_sizes.erase(the_model.subtree_sizes.begin()+2);
    the_model.node_names.erase(the_model.node_names.begin()+2);
    the_model.parent_ids.erase(the_model.parent_ids.begin()+2);
    the_model.branches.erase(the_model.branches.begin()+2);
    
    if (VERBOSE)
      cerr << the_model << endl;
    

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
  
    const size_t n_nodes = the_model.subtree_sizes.size();
    vector<vector<vector<double> > > all_interval_rates(n_nodes);
    vector<vector<double> > all_interval_lengths(n_nodes);
    rates_on_branch(the_model.triplet_rates,
                    all_paths[1][test_site - 1],
                    all_paths[1][test_site + 1],
                    all_interval_rates[1],
                    all_interval_lengths[1]);
    
    // (1) Upward pruning
    vector<vector<vector<double> > > all_p;
    all_p.resize(the_model.subtree_sizes.size());
    pruning(the_model.triplet_rates, the_model.subtree_sizes, test_site,
            all_paths, all_interval_rates, all_interval_lengths, all_p);
    
    // Downward sampling
    vector<Path> proposed_path_fs;
    size_t progress = 0;
    
    //while (fs_summary.size() < n_paths_to_sample) {
      if (VERBOSE && proposed_path_fs.size() * 10 % n_paths_to_sample == 0) {
        cerr << "FINISHED: " << progress << '%' << endl;
        progress += 10;
      }
      
      Path new_path, new_path_fs;
      
      // sample new root state
      const size_t root_id = 0; //all_paths[0] is empty
      const double root_p0 = root_post_prob0(test_site, all_paths,
                                             the_model.init_T, all_p);
      std::uniform_real_distribution<double> unif(0.0, 1.0);
         bool new_root_state = (unif(gen) > root_p0);
    
      new_path_fs.init_state = new_root_state;
      new_path_fs.tot_time = all_paths[1][test_site - 1].tot_time;
   
      // (2) Foward sample a path conditioned on given leaf state
    /*
      downward_sampling_branch_fs(all_interval_rates[1],
                                  all_interval_lengths[1],
                                  test_site,
                                  all_paths[1], all_p[1],
                                  gen, new_path_fs, max_iterations);
  
     */

    
      new_path.init_state = new_root_state;
      new_path.tot_time = all_paths[1][test_site - 1].tot_time;
      /*
      // (3) Downward sampling: pruning
      vector<Path> new_path_ds_all;
      downward_sampling(the_model.triplet_rates, the_model.subtree_sizes,
                        test_site, all_paths, the_model.init_T, all_p, gen,
                        new_path_ds_all);
      
      SummarySet current_summary_ds(all_paths[test_branch][test_site-1],
                                    new_path_ds_all[test_branch],
                                    all_paths[test_branch][test_site+1],
                                    n_hist_time_bins);
      ds_summary.push_back(current_summary_ds);
      
           SummarySet current_summary_fs(all_paths[test_branch][test_site-1],
                                    new_path_fs,
                                    all_paths[test_branch][test_site+1],
                                    n_hist_time_bins);
      fs_summary.push_back(current_summary_fs);
    }
    cerr << "FINISHED: " << progress << '%' << endl;
    
    // get the summaries of the summaries
    SummaryStatsFreq FS_report(fs_summary);
    SummaryStatsFreq DS_report(ds_summary);
    
    // write output
    std::ofstream out(outfile.c_str());
    out << "TEST_SITE_BRANCH" << '\t' << test_site << '\t'
                              << test_branch << endl;
    out << "N_PATHS"    << '\t' << n_paths_to_sample << endl;
    out << "IN_STATE" << '\t'
        << "PVAL_JUMPS" << '\t' << "PVAL_TIME" << '\t'
        << "HIST_JUMPS_FS" << '\t' << "HIST_JUMPS_DS" << '\t'
        << "HIST_TIME_FS" << '\t' << "HIST_TIME_DS" << endl;
    
    for(size_t k = 0; k < N_TRIPLETS; ++k) {
      double pval_time, pval_jumps;
      test_summary(FS_report, DS_report, pval_time, pval_jumps, k);
      out << std::bitset<3>(k).to_string() << '\t'
          << pval_jumps << '\t' << pval_time << '\t'
          << FS_report.print_jumps(k) << '\t'
          << DS_report.print_jumps(k) << '\t'
          << FS_report.print_time(k) << '\t'
          << DS_report.print_time(k) << endl;
    }
    
    }*/
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
