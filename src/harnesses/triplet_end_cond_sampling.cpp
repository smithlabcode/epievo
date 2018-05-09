/* Copyright (C) 2018 University of Southern California
 *                    Jianghan Qu and Andrew D Smith
 *
 * Author: Jianghan Qu, Andrew D. Smith and Xiaojing Ji
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
#include <numeric> /* iota, accumulate */
#include <algorithm>
#include <bitset>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "PhyloTreePreorder.hpp"
#include "Path.hpp"  /* related to Path */
#include "EpiEvoModel.hpp" /* model_param */
#include "StateSeq.hpp"
#include "SingleSampler.hpp"
#include "EndCondSampling.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::min;
using std::runtime_error;
using std::bitset;

static const double TIME_TOL = 1e-8;

bool
file_is_readable(const string &param_file) {
  std::ifstream in(param_file.c_str());
  return in.good();
}


static void
stats_from_jumps(const size_t is, const double tot_time,
                 const vector<double> &mid_jumps,
                 size_t &J, double &D) {

  Path l(get_left_bit(is), tot_time);
  Path r(get_right_bit(is), tot_time);
  Path m(get_mid_bit(is), tot_time, mid_jumps);

  PathContextStat pcs(l, m, r);
  J = pcs.jumps_in_context[is];
  D = pcs.time_in_context[is];
}

/* This function is copied from epievo_sim*/
static void
sample_jump_mid(const EpiEvoModel &the_model,
                const size_t is, const double tot_time,
                std::mt19937 &gen, vector<double> &jump_times,
                double &time_value) {

  // holding_rate = c_{ijk}*lambda_{ijk}
  const double holding_rate = the_model.triplet_rates[is];

  // sample a holding time = time until next state change
  std::exponential_distribution<double> exp_distr(holding_rate);
  const double holding_time = std::max(exp_distr(gen), TIME_TOL);
  // update the current time_value
  time_value += holding_time;

  // if the holding time ends before the total time interval, we can
  // make a change to the state sequence
  if (time_value < tot_time) {
    /* add the changed position and change time to the path */
    jump_times.push_back(time_value);
  }
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char **argv) {
  try {

    bool VERBOSE = false;
    bool SCALE = true;
    string outfile;

    size_t num_samples = 1000;
    size_t branch_to_sample = 1;

    string param_file;
    string tree_file;

    size_t rng_seed = std::numeric_limits<size_t>::max();

    ///////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]), "test sampling",
                           " <paths-file>");
    opt_parse.add_opt("param", 'p', "parameter file",
                      true, param_file);
    opt_parse.add_opt("tree", 't', "tree file in newick format",
                      true, tree_file);
    opt_parse.add_opt("branch", 'b', "branch to sample",
                      false, branch_to_sample);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);
    opt_parse.add_opt("number", 'n', "number of samples for each init state",
                      false, num_samples);
    opt_parse.add_opt("outfile", 'o', "outfile (sampling statistics) prefix",
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
    /*
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
     */

    if (!file_is_readable(param_file)) {
      cerr << "cannot read file: "<< param_file << endl
      << opt_parse.help_message() << endl;
      return  EXIT_SUCCESS;
    }

    ///////////////////////////////////////////////////////////////////////////
    /* standard mersenne_twister_engine seeded with rd() */
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    if (VERBOSE)
      cerr << "rng seed: " << rng_seed << endl;

    std::mt19937 gen(rng_seed);

    /* load parameters and tree */
    if (VERBOSE)
      cerr << "reading parameter file: " << param_file << endl;

    EpiEvoModel the_model;
    read_model(SCALE, param_file, tree_file, the_model);
    const double tot_time = the_model.branches[branch_to_sample];

    if (VERBOSE)
      cerr << the_model << endl;

    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    const size_t test_pos = 1; // test mid bit of triplet

    if (VERBOSE) {
      cerr << "total time: " << tot_time << endl;
      cerr << "number of samples for each init state: " << num_samples << endl;
    }

    cerr << "----- COMPARE forward / direct sampling ---------" << endl;
    vector<size_t> all_states_idx(8);
    std::iota(all_states_idx.begin()+1, all_states_idx.end(), 1);

    for(size_t i = 0; i < all_states_idx.size(); ++i) {
      const size_t i_bar = flip_mid_bit(i);
      const string init_triplet = bitset<3>(i).to_string();
      const string flip_triplet = bitset<3>(i_bar).to_string();

      if (VERBOSE)
          cerr <<  "---- TEST state: " << init_triplet << endl;

      const string outfile_summary = outfile + "." + init_triplet;
      std::ofstream outsummary(outfile_summary.c_str());

      vector<double> triplet_props;
      the_model.get_stationary_triplet_proportions(triplet_props);
      outsummary << "TOT_TIME" << '\t' << tot_time << '\t'
                 << 1.0 / the_model.triplet_rates[i] << endl;
      outsummary << "INIT_STATE" << '\t' << init_triplet << '\t'
                 << triplet_props[i] << endl;

      // (0) set up rates and helper variables
      vector<double> rates(2, 0.0);
      rates[0] = the_model.triplet_rates[i];
      rates[1] = the_model.triplet_rates[i_bar];
      vector<double> eigen_vals;
      vector<vector<double> > U;
      vector<vector<double> > Uinv;
      vector<vector<double> > PT;

      std::uniform_real_distribution<double> unif(0.0, 1.0);

      decompose(rates, eigen_vals, U, Uinv);
      continuous_time_trans_prob_mat(rates[0], rates[1], tot_time, PT);

      outsummary << "END_STATE_DISTR" << '\t'
                 << init_triplet << '\t' << PT[0][0] << '\t'
                 << flip_triplet << '\t' << PT[0][1] << endl;

      outsummary << "SAMPLE_IDX" << '\t' << "STATE_CHANGE" << '\t'
                 << "FS_SUCCEED" << '\t' << "NUM_FS_ATTEMPTS" << '\t'
                 << "J_FS" << '\t' << "D_FS" << '\t'
                 << "J_DS" << '\t' << "D_DS" << endl;

      //// statistics
      size_t num_success_samples = 0;
      size_t tot_fs_num_jumps = 0;
      double tot_fs_stay_time = 0;
      size_t tot_ds_num_jumps = 0;
      double tot_ds_stay_time = 0;

      for (size_t n = 0; n < num_samples; ++n) {
        // (1) sample an end_state according to transition rate
        const bool end_state_change = unif(gen) > PT[0][0];

        // (2) foward sampling
        const size_t MAX_ATTEMPTS = 1000;
        size_t num_sampled = 0;
        size_t curr_state = i;

        bool fs_reach_target = false;
        vector<double> fs_jump_times;

        while (!fs_reach_target && num_sampled < MAX_ATTEMPTS) {
          fs_jump_times.clear();

          double time_value = 0;
          while (time_value < tot_time) {
            sample_jump_mid(the_model, curr_state, tot_time, gen,
                            fs_jump_times, time_value);
            curr_state = flip_mid_bit(curr_state);
          }

          const bool sampled_change = curr_state == i;
          fs_reach_target = end_state_change == sampled_change;

          ++num_sampled;
        }

        size_t fs_num_jumps;
        double fs_stay_time;
        stats_from_jumps(i, tot_time, fs_jump_times, fs_num_jumps, fs_stay_time);


        // (3) direct sampling
        vector<double> ds_jump_times;
        size_t ds_num_jumps;
        double ds_stay_time;
        end_cond_sample(rates, eigen_vals, U, Uinv,
                        0, (size_t)(end_state_change),
                        tot_time, gen, ds_jump_times);
        stats_from_jumps(i, tot_time, ds_jump_times, ds_num_jumps, ds_stay_time);

        // (4) output line only when FS succeeds
        outsummary << n << '\t' << end_state_change << '\t'
                   << fs_reach_target << '\t' << num_sampled << '\t'
                   << fs_num_jumps << '\t' << fs_stay_time << '\t'
                   << ds_num_jumps << '\t' << ds_stay_time << endl;

        if (fs_reach_target) {
          ++num_success_samples;
          tot_fs_num_jumps += fs_num_jumps;
          tot_fs_stay_time += fs_stay_time;
          tot_ds_num_jumps += ds_num_jumps;
          tot_ds_stay_time += ds_stay_time;
        }
      }
      if (VERBOSE)
        cerr << "[NUM_SUCCESS_SAMPLES] " << num_success_samples
        << " [AVG_FS_D] " << tot_fs_num_jumps / num_success_samples
        << " [AVG_FS_D] " << tot_fs_stay_time / num_success_samples
        << " [AVG_FS_D] " << tot_ds_num_jumps / num_success_samples
        << " [AVG_FS_D] " << tot_ds_stay_time / num_success_samples
        << endl;
    }
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
