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
#include <algorithm>
#include <bitset>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "Path.hpp"  /* related to Path */
#include "EpiEvoModel.hpp" /* model_param */
#include "StateSeq.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::runtime_error;


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
                           " <paths-file>");
    opt_parse.add_opt("param", 'p', "parameter file",
                      true, param_file);
    opt_parse.add_opt("evo-time", 't', "time to simulate",
                      false, evo_time);
    opt_parse.add_opt("paths", 'p', "number of paths to sample",
                      false, n_paths_to_sample);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);
    opt_parse.add_opt("outfile", 'o', "outfile for posterior-updated paths",
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
    for (size_t i = 0; i < n_pairs; ++i) {

      // extract the left and right states
      bool left_state = false, right_state = false;
      get_bits_from_pair(i, left_state, right_state);

      // iterate over (both) start points for the mid state
      for (size_t j = 0; j < 2; ++j) {
        const bool mid_state = (i & 1ul);

        // simulate to obtain the desired number of paths
        vector<SummarySet> summary0, summary1;
        size_t iter = 0;
        while (iter++ < max_iterations &&
               summary0.size() < n_paths_to_sample &&
               summary1.size() < n_paths_to_sample) {

          // simulate a path starting at mid_state using forward simulation

          // check end state
          if (end_state && summary1.size() < n_paths_to_sample) {
            SummarySet current_summary;
            // obtain the summary stats; but only if still needed
            summary1.push_back(current_summary);
          }
          if (!end_state && summary0.size() < n_paths_to_sample) {
            SummarySet current_summary;
            // obtain the summary stats; but only if still needed
            summary0.push_back(current_summary);
          }
        }
        // now get the summaries of the summaries...?
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
