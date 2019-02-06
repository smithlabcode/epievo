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
#include <math.h>       /* floor */

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
using std::exp;


///////////////////////////////////////////////////////////////////////////
inline size_t
sequence2idx(const vector<char> &sequence) {
  assert(sequence.size() > 0);
  size_t idx = 0;
  size_t multiplier = 1;
  for (size_t i = sequence.size(); i > 0; i--) {
    idx += multiplier * (sequence[i-1] == true? 1 : 0);
    multiplier *= 2;
  }
  return idx;
}


///////////////////////////////////////////////////////////////////////////
//// "summary of pentets holding time"

struct PentetStat {
  PentetStat(const double _tot_time, const size_t n_bin);
  
  size_t time_to_bin_idx(const double time) const;
  void update_stat(const vector<size_t> starting_state,
                   const vector<double> jump_times);
  void update_stat(const char ll, const char l, const char m,
                   const char r, const char rr, const double time);
  string print_summary() const;
  
  double tot_time;
  vector<vector<size_t> > holding_time_count;
  vector<double> holding_time_sum;
  vector<size_t> pentet_count;
  double time_window;
  //vector<double> tw; // time window
  // histogram: (tw[0], tw[1]), ... (tw.back(), tot_time)
};

PentetStat::PentetStat(const double _tot_time, const size_t n_bin) {
  tot_time = _tot_time;
  
  holding_time_count = vector<vector<size_t> > (32, vector<size_t> (n_bin, 0));
  pentet_count = vector<size_t> (32, 0);
  holding_time_sum = vector<double> (32, 0.0);
  time_window = tot_time / n_bin;
}

size_t
PentetStat::time_to_bin_idx(const double time) const {
  assert(time > 0 && time < tot_time);
  //cerr << "time is: " << time << "time window is: " << time_window << endl;
  return static_cast<size_t>(floor( time / time_window));
}

void
PentetStat::update_stat(const vector<size_t> starting_state,
                        const vector<double> jump_times) {
  const size_t n_jumps = jump_times.size();
  assert(jump_times.size() == starting_state.size());
  
  // Note: last interval is not useful (overshoot)
  double last_jump_time = 0;
  for (size_t i = 0; i < n_jumps; i++) {
    const double time_interval = jump_times[i] - last_jump_time;
    holding_time_count[starting_state[i]][time_to_bin_idx(time_interval)] ++;
    pentet_count[starting_state[i]]++;
    holding_time_sum[starting_state[i]] += time_interval;
    
    last_jump_time = jump_times[i];
  }
}

void
PentetStat::update_stat(const char ll, const char l, const char m,
                        const char r, const char rr, const double time) {
  const vector<char> sequence = {ll, l, m, r, rr};
  const size_t pentet_idx = sequence2idx(sequence);
  holding_time_count[pentet_idx][time_to_bin_idx(time)] ++;
  pentet_count[pentet_idx]++;
  holding_time_sum[pentet_idx] += time;
}

string
PentetStat::print_summary() const {

  string str = "";
  for (size_t i = 0; i < holding_time_count.size(); i++) {
    const string pent_state = std::bitset<5>(i).to_string();
    str += pent_state + "," + std::to_string(pentet_count[i] > 0 ?
                                       holding_time_sum[i] / pentet_count[i] : 0);
    if (holding_time_count.size() > 0)
      str += "," + std::to_string(pentet_count[i] > 0 ?
           static_cast<double>(holding_time_count[i][0]) / pentet_count[i] : 0);
    for(size_t j = 1; j < holding_time_count[i].size(); j++) {
      str += "," + std::to_string(pentet_count[i] > 0 ?
           static_cast<double>(holding_time_count[i][j]) / pentet_count[i] : 0);
    }
    str += "\n";
  }
  return str;
}


///////////////////////////////////////////////////////////////////////////

static void
estimate_pentet_rate(const EpiEvoModel &model, const PentetStat pent,
                     vector<double> &estimated_rates,
                     vector<double> &my_guess) {
  for (size_t i = 0; i < pent.holding_time_sum.size(); i++) {
    if (pent.pentet_count[i] > 0)
      estimated_rates[i] = pent.pentet_count[i] / pent.holding_time_sum[i];

    const string pent_state = std::bitset<5>(i).to_string();
    const bool ll = (pent_state.at(0) ? true : false);
    const bool l = (pent_state.at(1) ? true : false);
    const bool m = (pent_state.at(2) ? true : false);
    const bool r = (pent_state.at(3) ? true : false);
    const bool rr = (pent_state.at(4) ? true : false);

    const size_t triplet_idx = triple2idx(l, m, r);
    my_guess[i] = model.triplet_rates[triplet_idx] *
                  exp(model.stationary_logbaseline[ll][m]) *
                  exp(model.stationary_logbaseline[ll][complement_state(m)]) *
                  exp(model.stationary_logbaseline[m][rr]) *
                  exp(model.stationary_logbaseline[complement_state(m)][rr]);
  }
}

///////////////////////////////////////////////////////////////////////////

/* This function does the sampling for an individual change in the
   state sequence
 */
static void
sample_jump(const EpiEvoModel &the_model, const double total_time,
            std::mt19937 &gen, TripletSampler &ts, vector<GlobalJump> &the_path,
            double &time_value, PentetStat &pent) {
  
  vector<char> seq;
  ts.get_sequence(seq);
  
  static const size_t n_triplets = 8;
  
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

    vector<double> triplet_prob(n_triplets, 0.0);
    for (size_t i = 0; i < n_triplets; ++i)
      triplet_prob[i] =
        triplet_counts[i]*the_model.triplet_rates[i]/holding_rate;

    std::discrete_distribution<size_t> multinom(triplet_prob.begin(),
                                                triplet_prob.end());
    const size_t context = multinom(gen);

    /* sample a change position having the relevant triplet; this
       changes the TripletSampler data structure to reflect a changed
       state at the position sampled */
    const size_t change_position = ts.random_mutate(context, gen);
    
    
    // update pentet stats
    if (change_position > 1 && change_position < seq.size() - 2)
      pent.update_stat(seq[change_position - 2], seq[change_position - 1],
                       seq[change_position], seq[change_position + 1 ],
                       seq[change_position + 2], holding_time);

    ts.get_sequence(seq);
    /* add the changed position and change time to the path */
    the_path.push_back(GlobalJump(time_value, change_position));
  }
}


///////////////////////////////////////////////////////////////////////////


int main(int argc, const char **argv) {
  try {

    static const size_t max_sample_count = 10000;

    bool VERBOSE = false;
    string outfile;

    size_t n_paths_to_sample = 1000;
    size_t n_hist_time_bins = 10;
      
    size_t n_sites = 5;
    double evo_time = 1.0;

    size_t rng_seed = std::numeric_limits<size_t>::max();

    ///////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]), "report pentet holding time",
                           "<param-file>");
    opt_parse.add_opt("sites", 's', "number of sites", false, n_sites);
    opt_parse.add_opt("time", 't', "duration of time interval", false, evo_time);
    opt_parse.add_opt("samples", 'n', "number of paths to sample",
                      false, n_paths_to_sample);
    opt_parse.add_opt("bins", 'b', "number of bins of holding time histogram",
                      false, n_hist_time_bins);
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
    //the_model.scale_triplet_rates();
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
    
    vector<char> root_seq;
    the_model.sample_state_sequence_init(n_sites, gen, root_seq);
    
    StateSeq s(root_seq);

    if (VERBOSE)
      cerr << "[SUMMARY:]" << endl
           << s.summary_string() << endl;


    PentetStat pentet_stat = PentetStat(evo_time, n_hist_time_bins);
    /* (4) SAMPLE THE PATHS GLOBALLY */
    for (size_t i = 0; i < n_paths_to_sample; ++i) {
      TripletSampler ts(s);

      vector<GlobalJump> the_path;
      double time_value = 0;

      while (time_value < evo_time)
        sample_jump(the_model, evo_time, gen, ts, the_path, time_value,
                    pentet_stat);
    }
    
    if (VERBOSE)
      cerr << "[Pentet holding time distribution:]" << endl;
    
    std::ofstream out(outfile.c_str());
    out << pentet_stat.print_summary() << endl;
    /*
    cerr << "Estimate pentet rates:" << endl;
    vector<double> estimated_pentet_rates = vector<double>(32, 0.0);
    vector<double> my_guesses = vector<double>(32, 0.0);
    estimate_pentet_rate(the_model, pentet_stat, estimated_pentet_rates,
                         my_guesses);
    for (size_t i = 0; i < estimated_pentet_rates.size(); ++i) {
      const string pent_state = std::bitset<5>(i).to_string();
      cerr << pent_state << ", " << estimated_pentet_rates[i]
      << ", " << my_guesses[i] << endl;
    }
    */
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
