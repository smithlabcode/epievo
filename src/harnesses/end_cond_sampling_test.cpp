/* Copyright (C) 2018 University of Southern California
 *                    Jianghan Qu and Andrew D Smith
 *
 * end_cond_sampling_test: This program is to test the methods for
 * end-conditioned sampling. We assume that the "sequence" only
 * considers two neighboring sites, and uses forward simulation as the
 * standard against which to test other methods. For each start-end
 * pair of states (binary) we want to compare a collection of summary
 * statistics obtained from sampling methods to those that would be
 * obtained by forward sampling.
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
#include <algorithm> /* count */
#include <bitset>
#include <gsl/gsl_randist.h> /* chi-squared test */

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "Path.hpp"  /* related to Path */
#include "EpiEvoModel.hpp" /* model_param */
#include "StateSeq.hpp"

#include <gsl/gsl_cdf.h>

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::runtime_error;


static double
evaluate_histogram_fit(const vector<double> &reference,
                       const vector<double> &to_evaluate) {

  assert(reference.size() == to_evaluate.size());

  double chi_squared_stat = 0.0;
  for (size_t i = 0; i < reference.size(); ++i) {
    chi_squared_stat +=
      (to_evaluate[i] - reference[i])*
      (to_evaluate[i] - reference[i])/reference[i];
  }

  const double degrees_of_freedom = reference.size() - 1;
  return gsl_cdf_chisq_P(chi_squared_stat, degrees_of_freedom);
}


struct SummarySet {
  SummarySet(const size_t J, const double D,
             const double mean_D) :
  num_jumps(J), total_stay_time(D), mean_stay_time(mean_D) {}

  SummarySet(const vector<double> jump_times);

  size_t num_jumps;
  double total_stay_time;
  double mean_stay_time;
};

SummarySet::SummarySet(const vector<double> jumps) {
  for (size_t i = 1; i < jumps.size(); ++i) {
    ++num_jumps;
    total_stay_time += (jumps[i] - jumps[i-1]);
  }
  --num_jumps; // last break point is not a jump
  mean_stay_time = total_stay_time / ((int) num_jumps / 2 + 1);
}

struct SummaryStatsFreq {
  SummaryStatsFreq() {}
  SummaryStatsFreq(const size_t n, const vector<size_t> J_freq,
                   const vector<size_t> &mean_D_freq) :
    num_samples(n), jumps_freq(J_freq), stay_time_freq(mean_D_freq) {}

  SummaryStatsFreq(const vector<SummarySet> &summary);

  size_t num_samples;
  size_t jumps_binsize;
  double time_binsize;
  vector<size_t> jumps_freq;
  vector<size_t> stay_time_freq; // average stay time
};

SummaryStatsFreq::SummaryStatsFreq(const vector<SummarySet> &summary) {
}



/* Forward sampling mid bit*/
static void
sample_jump_mid(const EpiEvoModel &the_model,
                const size_t is, const double tot_time,
                std::mt19937 &gen, vector<double> &jump_times,
                double &time_value) {

  static const size_t n_triplets = 8;

  // holding_rate = c_{ijk}*lambda_{ijk}
  const double holding_rate = the_model.triplet_rates[is];

  // sample a holding time = time until next state change
  std::exponential_distribution<double> exp_distr(holding_rate);
  const double holding_time =
    std::max(exp_distr(gen), std::numeric_limits<double>::min());

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
    string outfile;
    string param_file;

    size_t max_iterations = 1000000;
    size_t n_paths_to_sample = 1000;

    double evo_time = 1.0;

    size_t rng_seed = std::numeric_limits<size_t>::max();

    ///////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]), "test end-conditioned samplers",
                           "");
    opt_parse.add_opt("param", 'p', "parameter file",
                      true, param_file);
    opt_parse.add_opt("evo-time", 't', "time to simulate",
                      false, evo_time);
    opt_parse.add_opt("paths", 'p', "number of paths to sample",
                      false, n_paths_to_sample);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    ///////////////////////////////////////////////////////////////////////////

    if (VERBOSE)
      cerr << "time to simulate: " << evo_time << endl;

    /* standard mersenne_twister_engine seeded with rd() */
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    if (VERBOSE)
      cerr << "rng seed: " << rng_seed << endl;
    std::mt19937 gen(rng_seed);

    if (VERBOSE)
      cerr << "[READING PARAMETER FILE: " << param_file << "]" << endl;
    EpiEvoModel the_model;
    // ADS: the "false" below is to not scale the model
    read_model(param_file, the_model);

    if (VERBOSE)
      cerr << the_model << endl;

    // iterate over the possible contexts (left and right)
    const size_t n_pairs = 4;
    for (size_t i = 0; i < n_pairs; ++i) {
      // extract the left and right states
      bool left_state = false, right_state = false;
      get_bits_from_pair(i, left_state, right_state);

      // iterate over (both) start points for the mid state
      for (size_t j = 0; j < 2; ++j) {
        const bool mid_state = (i & 1ul);
        size_t triplet_idx = triple2idx(left_state, mid_state,
                                              right_state);
        // simulate to obtain the desired number of paths
        vector<SummarySet> summary0, summary1;
        size_t iter = 0;
        while (iter++ < max_iterations &&
               summary0.size() < n_paths_to_sample &&
               summary1.size() < n_paths_to_sample) {

          // simulate a path starting at mid_state using forward simulation
          vector<double> fs_jump_times;
          double time_value = 0;
          while (time_value < tot_time) {
            sample_jump_mid(the_model, triplet_idx, evo_time, gen,
                            fs_jump_times, time_value);
            triplet_idx = triple2idx(left_state, !mid_state, right_state);
          }
          const bool end_state = fs_jump_times.size() % 2 == 0 ?
                                 mid_state : !mid_state;

          // check end state
          if (end_state && summary1.size() < n_paths_to_sample) {
            SummarySet current_summary(fs_jump_times);
            // obtain the summary stats; but only if still needed
            summary1.push_back(current_summary);
          }
          if (!end_state && summary0.size() < n_paths_to_sample) {
            SummarySet current_summary(fs_jump_times);
            // obtain the summary stats; but only if still needed
            summary0.push_back(current_summary);
          }
        }
        // now get the summaries of the summaries...?
        SummaryStatsFreq FS_report0;
        SummaryStatsFreq FS_report1;
      }
    }

    // next do the same thing for the end-conditioned samplers

    // now do the tests on the summaries vs. the forward sim

  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
