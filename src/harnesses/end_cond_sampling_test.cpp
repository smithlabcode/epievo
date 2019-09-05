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
#include <sstream>
#include <fstream>
#include <algorithm> /* max_element */
#include <bitset>
#include <random>
#include <functional>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "Path.hpp"  /* related to Path */
#include "EpiEvoModel.hpp" /* model_param */
#include "StateSeq.hpp"
#include "EndCondSampling.hpp"
#include "ContinuousTimeMarkovModel.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::to_string;
using std::runtime_error;

typedef vector<vector<double> > two_by_two;

static double
count_zero_jumps(const size_t start_state, const vector<double> &jump_times) {
  const size_t n_jumps = jump_times.size();
  if (n_jumps % 2 == 0)
    return n_jumps/2;
  else
    return (start_state == 0) ? (n_jumps - 1)/2 + 1 : (n_jumps - 1)/2;
}

static double
get_time_in_zero(const size_t start_state, const double T,
                 const vector<double> &jump_times) {

  if (jump_times.size() == 0)
    return (start_state == 0) ? T : 0.0;

  double time_in_zero = 0;
  size_t a = start_state;
  double prev_time = 0.0;
  for (size_t i = 0; i < jump_times.size(); ++i) {
    if (a == 0)
      time_in_zero += jump_times[i] - prev_time;
    prev_time = jump_times[i];
    a = complement_state(a);
  }
  if (a == 0)
    time_in_zero += T - prev_time;

  return time_in_zero;
}

struct summary_info {

  summary_info(const string &mn,
               const bool s, const bool e, const double et) :
  method_name(mn), start_state(s), end_state(e), evo_time(et) {}
  
  string method_name;
  bool start_state;
  bool end_state;
  double evo_time;

  vector<double> jump_counts;
  vector<double> time_in_zero;
  vector<double> zero_jumps;
  vector<double> zero_duration;
  vector<double> time_in_one;
  vector<double> one_jumps;
  vector<double> one_duration;
  vector<double> proposal_prob;
  
  vector<double> all_jumps;
  static bool report_jumps;
  static string summary_suffix;
  static string jumps_suffix;
  static string filename_prefix;

  void update(const vector<double> &jump_times) {
    jump_counts.push_back(jump_times.size());
    time_in_zero.push_back(get_time_in_zero(start_state, evo_time, jump_times));
    time_in_one.push_back(evo_time - time_in_zero.back());
    zero_jumps.push_back(count_zero_jumps(start_state, jump_times));
    one_jumps.push_back(jump_counts.back() - zero_jumps.back());
    zero_duration.push_back(time_in_zero.back()/std::max(1.0, zero_jumps.back()));
    one_duration.push_back(time_in_one.back()/std::max(1.0, one_jumps.back()));
    if (report_jumps)
      all_jumps.insert(end(all_jumps), begin(jump_times), end(jump_times));
  }
  
  void update(const vector<double> &jump_times, const double prob) {
    update(jump_times);
    proposal_prob.push_back(prob);
  }

  void report() const;
  string tostring() const;
};
string summary_info::summary_suffix = "info";
string summary_info::jumps_suffix = "jumps";
string summary_info::filename_prefix = "";
bool summary_info::report_jumps = false;


template <typename T> double
get_mean(const T &x) {
  return accumulate(begin(x), end(x), 0.0)/x.size();
}


void
summary_info::report() const {
  if (report_jumps) {
    const string jumps_filename =
    filename_prefix + method_name + "." + jumps_suffix;
    std::ofstream j_out(jumps_filename);
    if (!j_out)
      throw runtime_error("cannot write to: " + jumps_filename);
    j_out << "start" << '\t'
    << "end" << '\t'
    << "time" << '\t'
    << "jumps" << endl;
    const bool end_state = (all_jumps.size() % 2 == 0) ? start_state : !start_state;
    j_out << start_state << '\t' << end_state << '\t' << evo_time;
    for (auto &&i : all_jumps)
      j_out << "\t" << i;
    j_out << endl;
  }
  
  const string info_filename =
  filename_prefix + method_name + "." + summary_suffix;
  std::ofstream out(info_filename);
  if (!out)
    throw runtime_error("cannot write to: " + info_filename);
  out << tostring() << endl;
}


string
summary_info::tostring() const {

  const double mean_zero_time = get_mean(time_in_zero);
  const double mean_one_time = get_mean(time_in_one);
  const double mean_zero_duration = get_mean(zero_duration);
  const double mean_one_duration = get_mean(one_duration);
  const double mean_zero_jumps = get_mean(zero_jumps);
  const double mean_one_jumps = get_mean(one_jumps);
  const double mean_proposal_prob = get_mean(proposal_prob);

  std::ostringstream oss;
  oss.precision(3);

  oss << method_name << '\t' << start_state << '\t' << end_state << '\t'
  << mean_zero_jumps << '\t' << mean_one_jumps << '\t'
  << mean_zero_time << '\t' << mean_one_time << '\t'
  << mean_zero_duration << '\t' << mean_one_duration << '\t'
  << mean_proposal_prob << endl;
  
  return oss.str();
}



////////////////////////////////////////////////////////////////////////////////

template <typename T>
static void
append_to_file(const string &outfile, const T &x) {
  std::ofstream out(outfile, std::ofstream::app);
  out << x;
}


static string
expected_stat_str(const bool start_state, const bool end_state,
                  const two_by_two &J0, const two_by_two &J1,
                  const two_by_two &D0, const two_by_two &D1) {
  std::ostringstream oss;
  oss.precision(3);
  
  oss << start_state << '\t' << end_state << '\t'
  << J0[start_state][end_state] << '\t' << J1[start_state][end_state] << '\t'
  << D0[start_state][end_state] << '\t' << D1[start_state][end_state] << '\t'
  << "\\\t\\\t\\\t" << endl;
  
  return oss.str();
}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////




int main(int argc, const char **argv) {
  try {

    bool VERBOSE = false;
    string statfile;

    size_t n_samples = 1000;

    double rate0 = 1.5;
    double rate1 = 0.5;
    double evo_time = 1.0;

    size_t rng_seed = std::numeric_limits<size_t>::max();

    ///////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]), "test end-conditioned samplers",
                           "<output>");
    opt_parse.add_opt("rate0", '\0', "transition rate from state 0", false,
                      rate0);
    opt_parse.add_opt("rate1", '\0', "transition rate from state 1", false,
                      rate1);
    opt_parse.add_opt("time", 't', "duration of time interval", false, evo_time);
    opt_parse.add_opt("jumps", 'j', "write the jump times",
                      false, summary_info::report_jumps);
    opt_parse.add_opt("n-samples", 'n', "number of samples to simulate", false,
                      n_samples);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);
    opt_parse.add_opt("statfile", 'S', "summary file", false, statfile);

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
    
    const string outfile(leftover_args.front());
    ///////////////////////////////////////////////////////////////////////////

    /* INITIALIZING PARAMETERS */
    const TwoStateCTMarkovModel ctmm(rate0, rate1);
    if (VERBOSE)
      cerr << ctmm << endl;

    /* standard mersenne_twister_engine seeded with rd() */
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    if (VERBOSE)
      cerr << "RNG seed: " << rng_seed << endl << endl;
    std::mt19937 gen(rng_seed);

    
    const string statfile_header = "method\tstart\tend\tJ0\tJ1\tD0\tD1\ttau0\ttau1\tprob";
    if (VERBOSE)
      cerr << statfile_header << endl;
    if (!statfile.empty()) {
      std::ofstream outstat(statfile.c_str());
      outstat << statfile_header << endl;
    }

    two_by_two expected_J0, expected_J1, expected_D0, expected_D1;
    expectation_J(rate0, rate1, evo_time, expected_J0, expected_J1);
    expectation_D(rate0, rate1, evo_time, expected_D0, expected_D1);

    vector<double> jump_times;

    for (size_t start_state = 0; start_state <= 1; start_state++) {
      for (size_t end_state = 0; end_state <= 1; end_state++) {
        const string expstat = "*Expected     \t" +
        expected_stat_str(start_state, end_state, expected_J0, expected_J1,
                          expected_D0, expected_D1);
        if (VERBOSE)
          cerr << expstat;
        if(!statfile.empty())
          append_to_file(statfile, expstat);

        summary_info si_direct("direct", start_state, end_state, evo_time);
        summary_info si_forward("forward", start_state, end_state, evo_time);
        summary_info si_nielsen("nielsen", start_state, end_state, evo_time);
        summary_info si_unif("unif", start_state, end_state, evo_time);
        summary_info si_pois("pois", start_state, end_state, evo_time);
        
        /* Direct sampling */
        for (size_t i = 0; i < n_samples; i++) {
          jump_times.clear();
          end_cond_sample_direct(ctmm, start_state, end_state, evo_time, gen,
                                 jump_times);
        }
        si_direct.report();
        if (VERBOSE)
          cerr << si_direct.tostring();

        /* Forward + rejection sampling */
        for (size_t i = 0; i < n_samples; i++) {
          jump_times.clear();
          end_cond_sample_forward_rejection(ctmm, start_state, end_state,
                                            evo_time, gen, jump_times);
        }
        si_forward.report();
        if (VERBOSE)
          cerr << si_forward.tostring();

        /* Forward + rejection (Nielsen) sampling */
        for (size_t i = 0; i < n_samples; i++) {
          jump_times.clear();
          end_cond_sampling_Nielsen(ctmm, start_state, end_state,
                                    evo_time, gen, jump_times);
        }
        si_nielsen.report();
        if (VERBOSE)
          cerr << si_nielsen.tostring();
        
        /* Uniformization sampling */
        for (size_t i = 0; i < n_samples; i++) {
          jump_times.clear();
          end_cond_sample_unif(ctmm, start_state, end_state, evo_time, gen,
                               jump_times);
        }
        si_unif.report();
        if (VERBOSE)
          cerr << si_unif.tostring();

        /* Poisson sampling */
        for (size_t i = 0; i < n_samples; i++) {
          jump_times.clear();
          end_cond_sample_Poisson(ctmm, start_state, end_state, evo_time, gen,
                                  jump_times);
        }
        si_pois.report();
        if (VERBOSE)
          cerr << si_pois.tostring();
      }
    }
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
