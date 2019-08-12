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
#include <algorithm> /* max_element */
#include <bitset>
#include <random>
#include <functional>

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
#include "ContinuousTimeMarkovModel.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::runtime_error;


static const double MINWAIT = 1e-8;


static double
integral_J_ij(const double lambda_i, const double lambda_j, const double T) {
  return (lambda_i == lambda_j ?
          T*exp(lambda_i*T) :
          (exp(lambda_i*T) - exp(lambda_j*T))/(lambda_i - lambda_j));
}


static void
expected_time_in_state(const double rate0,
                       const double rate1,
                       const double T,
                       vector<double> &ED) {

  // here corresponding to the "a", "b" and "j"
  static const double n_triplets = 8;

  ED = vector<double>(n_triplets, 0.0);

  const vector<double> rates = {rate0, rate1};

  vector<double> lambda;
  vector<vector<double> > U;
  vector<vector<double> > Uinv;
  decompose(rates, lambda, U, Uinv);

  vector<vector<double> > P;
  continuous_time_trans_prob_mat(rate0, rate1, T, P);

  // compute integral_J_ij
  vector<vector<double>> J = vector<vector<double>>(2, vector<double>(2, 0.0));
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      J[i][j] = integral_J_ij(lambda[i], lambda[j], T);

      
  for (size_t start_state = 0; start_state < 2; ++start_state) {
    for (size_t end_state = 0; end_state < 2; ++end_state) {
      for (size_t c = 0; c < 2; ++c) {
        ED[triple2idx(start_state, end_state, c)] =
          (U[start_state][0]*Uinv[0][c]*(U[c][0]*Uinv[0][end_state]*J[0][0] +
                                         U[c][1]*Uinv[1][end_state]*J[0][1]) +
           U[start_state][1]*Uinv[1][c]*(U[c][0]*Uinv[0][end_state]*J[1][0] +
                                         U[c][1]*Uinv[1][end_state]*J[1][1]))/
          P[start_state][end_state];
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



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
  gsl_histogram_set_ranges_uniform(h_time, 0, tot_time+MINWAIT);

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
  }
  else {
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


////////////////////////////////////////////////////////////////////////////////

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

static double
mean_seg_pair_length(const double T,
                     const vector<double> &jump_times) {

  if (jump_times.size() % 2 == 0)
    return 0.0;

  double total_length = 0.0;
  size_t count = 0;
  double prev_time = 0.0;
  for (size_t i = 1; i < jump_times.size(); i += 2) {
    total_length += jump_times[i] - prev_time;
    prev_time = jump_times[i];
    ++count;
  }
  total_length += (T - prev_time);
  ++count;

  return total_length/count;
}


static double
count_zero_segments(const size_t start_state, const vector<double> &jump_times) {
  const size_t n_segs = jump_times.size() + 1;
  if (n_segs % 2 == 0)
    return n_segs/2;
  else
    return (start_state == 0) ? (n_segs - 1)/2 + 1 : (n_segs - 1)/2;
}

int main(int argc, const char **argv) {
  try {

    static const size_t max_sample_count = 10000;

    bool VERBOSE = false;
    string outfile;

    size_t n_paths_to_sample = 1000;
    size_t n_hist_time_bins = 10;

    double rate0 = 1.5;
    double rate1 = 0.5;
    double evo_time = 1.0;

    size_t rng_seed = std::numeric_limits<size_t>::max();

    ///////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]), "test end-conditioned samplers",
                           "");
    opt_parse.add_opt("rate0", '\0', "holding rate for state 0", true, rate0);
    opt_parse.add_opt("rate1", '\0', "holding rate for state 1", true, rate1);
    opt_parse.add_opt("time", 't', "duration of time interval", true, evo_time);
    opt_parse.add_opt("samples", 'n', "number of paths to sample",
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

    if (VERBOSE)
      cerr << "rate0: " << rate0 << endl
           << "rate1: " << rate1 << endl
           << "time to simulate: " << evo_time << endl
           << "paths to simulate: " << n_paths_to_sample << endl;

    /* standard mersenne_twister_engine seeded with rd() */
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    if (VERBOSE)
      cerr << "rng seed: " << rng_seed << endl;
    std::mt19937 gen(rng_seed);

    const CTMarkovModel ctmm(rate0, rate1);

    cout << ctmm << endl;

    vector<double> ED;
    expected_time_in_state(rate0, rate1, evo_time, ED);

    // direct sampling
    cout << "direct sampling:" << endl;
    for (size_t start_state = 0; start_state < 2; ++start_state) {
      for (size_t end_state = 0; end_state < 2; ++end_state) {

        vector<double> seg_counts;
        vector<double> time_in_zero;
        vector<double> zero_segs;
        vector<double> zero_duration;
        vector<double> time_in_one;
        vector<double> one_segs;
        vector<double> one_duration;
        vector<double> proposal_prob;
        for (size_t k = 0; k < n_paths_to_sample; ++k) {
          vector<double> jump_times;
          vector<mixJump> mjumps;
          //end_cond_sample_direct(ctmm, start_state, end_state, evo_time,
          //                       gen, jump_times);
          end_cond_sample_unif(ctmm, start_state, end_state, evo_time,
                               gen, jump_times, mjumps);
          seg_counts.push_back(jump_times.size() + 1);
          time_in_zero.push_back(get_time_in_zero(start_state, evo_time, jump_times));
          time_in_one.push_back(evo_time - time_in_zero.back());
          zero_segs.push_back(count_zero_segments(start_state, jump_times));
          one_segs.push_back(seg_counts.back() - zero_segs.back());
          // if (zero_segs.back() > 0)
          //   zero_duration.push_back(time_in_zero.back()/zero_segs.back());
          // if (one_segs.back() > 0)
          //   one_duration.push_back(time_in_one.back()/one_segs.back());
          zero_duration.push_back(time_in_zero.back()/std::max(1.0, zero_segs.back()));
          one_duration.push_back(time_in_one.back()/std::max(1.0, one_segs.back()));
          proposal_prob.push_back(end_cond_sample_prob(ctmm, jump_times,
                                                       start_state, end_state,
                                                       0, evo_time,
                                                       0, jump_times.size()));
        }
        const double mean_segs =
          accumulate(begin(seg_counts), end(seg_counts), 0.0)/seg_counts.size();

        const double mean_zero_time =
          accumulate(begin(time_in_zero), end(time_in_zero), 0.0)/time_in_zero.size();
        const double mean_one_time =
          accumulate(begin(time_in_one), end(time_in_one), 0.0)/time_in_one.size();

        const double mean_zero_duration =
          accumulate(begin(zero_duration), end(zero_duration), 0.0)/zero_duration.size();
        const double mean_one_duration =
          accumulate(begin(one_duration), end(one_duration), 0.0)/one_duration.size();

        const double mean_zero_segs =
          accumulate(begin(zero_segs), end(zero_segs), 0.0)/zero_segs.size();
        const double mean_one_segs =
          accumulate(begin(one_segs), end(one_segs), 0.0)/one_segs.size();

        const double mean_proposal_prob =
        accumulate(begin(proposal_prob), end(proposal_prob), 0.0)/proposal_prob.size();
        
        cout << "X(0)=" << start_state << '\t'
             << "X(T)=" << end_state << '\t'
             << "Segs=" << mean_segs << '\t'
             << "T(0)=" << mean_zero_time << '\t'
             << "S(0)=" << mean_zero_segs << '\t'
             << "D(0)=" << mean_zero_duration << '\t'
             << "T(1)=" << mean_one_time << '\t'
             << "S(1)=" << mean_one_segs << '\t'
             << "D(1)=" << mean_one_duration << '\t'
             << "Proposal prob=" << mean_proposal_prob << endl;
      }
    }

    cout << "naive forward sampling with rejection:" << endl;
    // naive forward sampling with rejection
    for (size_t start_state = 0; start_state < 2; ++start_state) {
      for (size_t end_state = 0; end_state < 2; ++end_state) {

        vector<double> seg_counts;
        vector<double> time_in_zero;
        vector<double> zero_segs;
        vector<double> zero_duration;
        vector<double> time_in_one;
        vector<double> one_segs;
        vector<double> one_duration;
        vector<double> proposal_prob;
        vector<double> mspl;
        for (size_t k = 0; k < n_paths_to_sample; ++k) {
          vector<double> jump_times;
          assert(end_cond_sample_forward_rejection(max_sample_count, ctmm,
                                                   start_state, end_state, evo_time,
                                                   gen, jump_times));
          seg_counts.push_back(jump_times.size() + 1);
          time_in_zero.push_back(get_time_in_zero(start_state, evo_time, jump_times));
          time_in_one.push_back(evo_time - time_in_zero.back());
          zero_segs.push_back(count_zero_segments(start_state, jump_times));
          one_segs.push_back(seg_counts.back() - zero_segs.back());
          // if (zero_segs.back() > 0)
          //   zero_duration.push_back(time_in_zero.back()/zero_segs.back());
          // if (one_segs.back() > 0)
          //   one_duration.push_back(time_in_one.back()/one_segs.back());
          zero_duration.push_back(time_in_zero.back()/std::max(1.0, zero_segs.back()));
          one_duration.push_back(time_in_one.back()/std::max(1.0, one_segs.back()));
          if (jump_times.size() % 2 == 1)
            mspl.push_back(mean_seg_pair_length(evo_time, jump_times));
          proposal_prob.push_back(forward_sample_prob(ctmm, jump_times,
                                                      start_state, end_state,
                                                      0, evo_time,
                                                      0, jump_times.size()));
        }
        const double mean_segs =
          accumulate(begin(seg_counts), end(seg_counts), 0.0)/seg_counts.size();

        const double mean_zero_time =
          accumulate(begin(time_in_zero), end(time_in_zero), 0.0)/time_in_zero.size();
        const double mean_one_time =
          accumulate(begin(time_in_one), end(time_in_one), 0.0)/time_in_one.size();

        const double mean_zero_duration =
          accumulate(begin(zero_duration), end(zero_duration), 0.0)/zero_duration.size();
        const double mean_one_duration =
          accumulate(begin(one_duration), end(one_duration), 0.0)/one_duration.size();

        const double mean_zero_segs =
          accumulate(begin(zero_segs), end(zero_segs), 0.0)/zero_segs.size();
        const double mean_one_segs =
          accumulate(begin(one_segs), end(one_segs), 0.0)/one_segs.size();

        const double mean_mspl =
          accumulate(begin(mspl), end(mspl), 0.0)/mspl.size();

        const double mean_proposal_prob =
        accumulate(begin(proposal_prob), end(proposal_prob), 0.0)/proposal_prob.size();
        
        cout << "X(0)=" << start_state << '\t'
             << "X(T)=" << end_state << '\t'
             << "Segs=" << mean_segs << '\t'
             << "J=" << mean_segs - 1 << '\t'
             << "T(0)=" << mean_zero_time << '\t'
             << "S(0)=" << mean_zero_segs << '\t'
             << "D(0)=" << mean_zero_duration << '\t'
             << "T(1)=" << mean_one_time << '\t'
             << "S(1)=" << mean_one_segs << '\t'
             << "D(1)=" << mean_one_duration << '\t'
             << "MSP=" << mean_mspl << '\t'
             << "E[T(0)]=" << ED[triple2idx(start_state, end_state, 0ul)] << '\t'
             << "Proposal prob=" << mean_proposal_prob << endl;

        typedef std::poisson_distribution<size_t> pois_distr;
        auto pois = bind(pois_distr(mean_zero_duration + mean_one_duration), ref(gen));
        vector<double> jump_counts;
        for (size_t k = 0; k < n_paths_to_sample; ++k) {
          jump_counts.push_back(pois());
        }
        cout << accumulate(begin(jump_counts),
                           end(jump_counts), 0.0)/jump_counts.size() << endl;

        auto pois2 = bind(pois_distr(evo_time*(1.0/rate0 + 1.0/rate1)), ref(gen));
        jump_counts.clear();
        for (size_t k = 0; k < n_paths_to_sample; ++k) {
          jump_counts.push_back(pois2());
        }
        cout << accumulate(begin(jump_counts),
                           end(jump_counts), 0.0)/jump_counts.size() << endl;
      }
    }
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
