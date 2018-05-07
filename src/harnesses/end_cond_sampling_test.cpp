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
#include <algorithm> /* max_element */
#include <bitset>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_randist.h> /* chi-squared test */
#include <gsl/gsl_cdf.h>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "Path.hpp"  /* related to Path */
#include "EpiEvoModel.hpp" /* model_param */
#include "StateSeq.hpp"
#include "EndCondSampling.hpp"

/* test TripletSampler
#include "TripletSampler.hpp"
#include "GlobalJump.hpp"
*/

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::runtime_error;


static double
evaluate_fit(const vector<double> &reference,
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


///////////////////////////////////////////////////////////////////////////

struct SummarySet {
  SummarySet(const vector<double> &jumps, const double tot_time,
             const size_t n_bins);

  size_t num_jumps;
  double total_stay_time;
  gsl_histogram * h_time;
};

SummarySet::SummarySet(const vector<double> &jumps, const double tot_time,
                       const size_t n_bins) {
  // initialize histogram
  h_time = gsl_histogram_alloc(n_bins);
  gsl_histogram_set_ranges_uniform(h_time, 0, tot_time+1);
  
  num_jumps = 0;   // number of jump-backs in init state
  total_stay_time = 0;
  
  // retrieve N-1 intervals/jumps in init state
  double elapsed_time = 0;
  for (size_t i = 1; i < jumps.size(); i += 2) {
    ++num_jumps;
    total_stay_time += (jumps[i-1] - elapsed_time);
    gsl_histogram_increment(h_time, jumps[i-1] - elapsed_time);
    elapsed_time = jumps[i];
  }
  
  // last interval
  if (jumps.size() % 2 == 0) {
    total_stay_time += (tot_time - elapsed_time);
    gsl_histogram_increment(h_time, tot_time - elapsed_time);
  } else {
    total_stay_time += (jumps.back() - elapsed_time);
    gsl_histogram_increment(h_time, jumps.back() - elapsed_time);
  }
}

///////////////////////////////////////////////////////////////////////////
//// "summary of summaries"

struct SummaryStatsFreq {
  SummaryStatsFreq(const vector<SummarySet> &summary);
  string print_time() const;
  string print_jumps() const;

  size_t num_samples;
  vector<double> time_freq; // proportion, not count
  vector<double> jumps_freq; // proportion, not count
};

SummaryStatsFreq::SummaryStatsFreq(const vector<SummarySet> &summary) {
  num_samples = summary.size();
  
  if (num_samples > 0) {
    // initialize time_freq
    const size_t nbins_time = summary.back().h_time->n;
    time_freq.resize(nbins_time, 0);
    
    // initialize jumps_freq
    auto it = max_element(summary.begin(), summary.end(),
                          [] (SummarySet const& s1, SummarySet const& s2) {
                            return s1.num_jumps < s2.num_jumps;
                          });
    
    jumps_freq.resize(it->num_jumps+1, 0);
    
    // merge summaries
    for (size_t i = 0; i < num_samples; i++) {
      for (size_t j = 0; j < nbins_time; j++) {
        time_freq[j] += summary[i].h_time->bin[j];
      }
      jumps_freq[summary[i].num_jumps]++;
    }
    
    // normalize to frequency
    const double sum_time = accumulate(time_freq.begin(), time_freq.end(), 0.0);
    const double sum_jumps = accumulate(jumps_freq.begin(), jumps_freq.end(),
                                        0.0);
    
    for (size_t i = 0; i < time_freq.size(); i++) {
      time_freq[i] /= sum_time;
    }
    for (size_t i = 0; i < jumps_freq.size(); i++) {
      jumps_freq[i] /= sum_jumps;
    }
  }
}

string
SummaryStatsFreq::print_time() const {
  string str;
  if (time_freq.size() > 0)
    str = std::to_string(time_freq[0]);
  for(size_t i = 1; i < time_freq.size(); i++) {
    str = str + "," + std::to_string(time_freq[i]);
  }
  return str;
}


string
SummaryStatsFreq::print_jumps() const {
  string str;
  if (jumps_freq.size() > 0)
    str = std::to_string(jumps_freq[0]);
  for(size_t i = 1; i < jumps_freq.size(); i++) {
    str = str + "," + std::to_string(jumps_freq[i]);
  }
  return str;
}

static void
test_summary(SummaryStatsFreq &a, SummaryStatsFreq &b,
             double &pval_time, double &pval_jumps) {
  // rule out zero expected values
  vector<double> time_exp, time_obs, jump_exp, jump_obs;
  for(size_t i = 0; i < std::min(a.time_freq.size(), b.time_freq.size()); i++) {
    if (a.time_freq[i] > 0) {
      time_exp.push_back(a.time_freq[i]);
      time_obs.push_back(b.time_freq[i]);
    }
  }
  pval_time = evaluate_fit(time_exp, time_obs);
  
  for(size_t i = 0; i < std::min(a.jumps_freq.size(), b.jumps_freq.size()); i++) {
    if (a.jumps_freq[i] > 0) {
      jump_exp.push_back(a.jumps_freq[i]);
      jump_obs.push_back(b.jumps_freq[i]);
    }
  }
  pval_jumps = evaluate_fit(jump_exp, jump_obs);
}


///////////////////////////////////////////////////////////////////////////


/* Forward sampling mid bit*/
static void
sample_jump_mid(const EpiEvoModel &the_model,
                const size_t is, const double tot_time,
                std::mt19937 &gen, vector<double> &jump_times,
                double &time_value) {
  
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

/* Forward sampling mid bit using TripletSampler
static void
sample_jump_mid(const EpiEvoModel &the_model, const double total_time,
                std::mt19937 &gen, TripletSampler &ts, vector<GlobalJump> &the_path,
                double &time_value) {
  
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
                                       std::numeric_limits<double>::min());
  time_value += holding_time;
 
  if (time_value < total_time) {
 
    size_t context;
    for(size_t i=0; i < triplet_counts.size(); i++) {
      if (triplet_counts[i] > 0)
        context = i;
    }
    
    const size_t change_position = ts.random_mutate(context, gen);
    the_path.push_back(GlobalJump(time_value, change_position));
  }
}
*/

////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char **argv) {
  try {

    bool VERBOSE = false;
    string outfile;
    string param_file;

    size_t max_iterations = 1000000;
    size_t n_paths_to_sample = 1000;
    size_t n_hist_time_bins = 10;
    
    double evo_time = 1.0;

    size_t rng_seed = std::numeric_limits<size_t>::max();

    ///////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]), "test end-conditioned samplers",
                           "");
    opt_parse.add_opt("param", 'p', "parameter file",
                      true, param_file);
    opt_parse.add_opt("evo-time", 't', "time to simulate",
                      false, evo_time);
    opt_parse.add_opt("paths", 'n', "number of paths to sample",
                      false, n_paths_to_sample);
    opt_parse.add_opt("bins", 'b', "number of bins of holding time histogram",
                      false, n_hist_time_bins);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);
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
    ///////////////////////////////////////////////////////////////////////////

    if (VERBOSE) {
      cerr << "time to simulate: " << evo_time << endl;
      cerr << "paths to simulate: " << n_paths_to_sample << endl;
    }

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

    std::ofstream out(outfile.c_str());
    out << "TOTAL_TIME" << '\t' << evo_time << endl;
    out << "N_PATHS"    << '\t' << n_paths_to_sample << endl;
    out << "INIT_STATE" << '\t' << "END_STATE" << '\t'
        << "PVAL_JUMPS" << '\t' << "PVAL_TIME" << '\t'
        << "HIST_JUMPS_FS" << '\t' << "HIST_JUMPS_DS" << '\t'
        << "HIST_TIME_FS" << '\t' << "HIST_TIME_DS" << endl;
    
    
    
    // iterate over the possible contexts (left and right)
    const size_t n_pairs = 4;
    for (size_t i = 0; i < n_pairs; ++i) {
      // extract the left and right states
      bool left_state = false, right_state = false;
      get_bits_from_pair(i, left_state, right_state);

      // iterate over (both) start points for the mid state
      for (size_t j = 0; j < 2; ++j) {
        const bool mid_state = (j & 1ul);
        size_t triplet_idx = triple2idx(left_state, mid_state,
                                        right_state);
        
        if (VERBOSE)
          cerr <<  "FORWARD SAMPLING INIT STATE: "
               << std::bitset<3>(triplet_idx).to_string() << endl;

        // 1. forward samplers
        vector<SummarySet> fs_summary0, fs_summary1;
        size_t iter = 0;
        while (iter++ < max_iterations &&
               fs_summary0.size() < n_paths_to_sample &&
               fs_summary1.size() < n_paths_to_sample) {

          // simulate a path starting at mid_state using forward simulation
          vector<double> fs_jump_times;
          double time_value = 0;
          size_t curr_state = triplet_idx;
          
          /* test TripletSampler
          vector<char> seq = {left_state, mid_state, right_state};
          TripletSampler ts(seq);
          vector<GlobalJump> the_path;
          */
          while (time_value < evo_time) {
            sample_jump_mid(the_model, curr_state, evo_time, gen,
                            fs_jump_times, time_value);
            curr_state = flip_mid_bit(curr_state);
            /* test TripletSampler
            sample_jump_mid(the_model, evo_time, gen, ts,
                            the_path, time_value);
            */
            
          }
          /* testTripletSampler
          ts.get_sequence(seq);
          const bool end_state = seq[1];
          for(size_t k = 0; k < the_path.size(); k++)
            fs_jump_times.push_back(the_path[k].timepoint);
          */
          
          const bool end_state = fs_jump_times.size() % 2 == 0 ?
                                 mid_state : !mid_state;

          // check end state
          if (end_state && fs_summary1.size() < n_paths_to_sample) {
            SummarySet current_summary(fs_jump_times, evo_time,
                                       n_hist_time_bins);
            fs_summary1.push_back(current_summary);
          }
          if (!end_state && fs_summary0.size() < n_paths_to_sample) {
            SummarySet current_summary(fs_jump_times, evo_time,
                                       n_hist_time_bins);
            fs_summary0.push_back(current_summary);
          }
        }
        // now get the summaries of the summaries...?
        SummaryStatsFreq FS_report0(fs_summary0);
        SummaryStatsFreq FS_report1(fs_summary1);
        
        // 2. end-conditioned samplers
        cerr <<  "DIRECT SAMPLING INIT STATE: "
             << std::bitset<3>(triplet_idx).to_string() << endl;
        
        vector<SummarySet> ds_summary0, ds_summary1;
        
        // 2.0 set up rates and helper variables
        vector<double> rates(2, 0.0);
        rates[0] = the_model.triplet_rates[triplet_idx];
        rates[1] = the_model.triplet_rates[flip_mid_bit(triplet_idx)];
        vector<double> eigen_vals;
        vector<vector<double> > U;
        vector<vector<double> > Uinv;
        vector<vector<double> > PT;

        decompose(rates, eigen_vals, U, Uinv);
        trans_prob_mat(rates[0], rates[1], evo_time, PT);
        
        while (ds_summary0.size() < n_paths_to_sample) {
          vector<double> ds_jump_times0, ds_jump_times1;
          
          end_cond_sample(rates, eigen_vals, U, Uinv,
                          0, mid_state, evo_time, gen, ds_jump_times0);
          SummarySet current_summary0(ds_jump_times0, evo_time,
                                     n_hist_time_bins);
          ds_summary0.push_back(current_summary0);
          
          end_cond_sample(rates, eigen_vals, U, Uinv,
                          0, !mid_state, evo_time, gen, ds_jump_times1);
          SummarySet current_summary1(ds_jump_times1, evo_time,
                                      n_hist_time_bins);
          ds_summary1.push_back(current_summary1);
        }

        SummaryStatsFreq DS_report0(ds_summary0);
        SummaryStatsFreq DS_report1(ds_summary1);
        
        // 3. statistical testing and output
        double pval_time0, pval_jumps0, pval_time1, pval_jumps1;
        test_summary(FS_report0, DS_report0, pval_time0, pval_jumps0);
        test_summary(FS_report1, DS_report1, pval_time1, pval_jumps1);
        out << std::bitset<3>(triplet_idx).to_string() << '\t'
        << std::bitset<3>(triplet_idx).to_string() << '\t'
        << pval_jumps0 << '\t' << pval_time0 << '\t'
        << FS_report0.print_jumps() << '\t'
        << DS_report0.print_jumps() << '\t'
        << FS_report0.print_time() << '\t'
        << DS_report0.print_time() << endl;
        
        out << std::bitset<3>(triplet_idx).to_string() << '\t'
        << std::bitset<3>(flip_mid_bit(triplet_idx)).to_string() << '\t'
        << pval_jumps1 << '\t' << pval_time1 << '\t'
        << FS_report1.print_jumps() << '\t'
        << DS_report1.print_jumps() << '\t'
        << FS_report1.print_time() << '\t'
        << DS_report1.print_time() << endl;
      }
    }
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
