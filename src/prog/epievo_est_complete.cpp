/* Copyright (C) 2018 University of Southern California
 *                    Jianghan Qu and Andrew D Smith
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
#include <random>
#include <math.h>   /* exp, sqrt, pow fabs*/
#include <numeric>  /* std::inner_product, accumulate*/
#include <functional>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "PhyloTreePreorder.hpp"
#include "path.hpp"  /* related to Path */
#include "TwoStateSeq.hpp"
#include "param.hpp" /* model_param */
#include "TripletPattern.hpp"
#include "jump.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;

void
add_sufficient_statistics(const Path &left, const Path &mid, const Path &right,
                          vector<double> &J, vector<double> &D) {

  assert(left.tot_time == mid.tot_time && mid.tot_time == right.tot_time);

  size_t triplet = triple2idx(left.init_state, mid.init_state, right.init_state);
  double prev_time = 0.0;
  size_t i = 0, j = 0, k = 0;
  while (i < left.jumps.size() && j < mid.jumps.size() && k < right.jumps.size())
    if (left.jumps[i] < min(mid.jumps[j], right.jumps[k])) {
      D[triplet] += left.jumps[i] - prev_time;
      J[triplet] += 1.0;
      prev_time = left.jumps[i];
      triplet = flip_left_bit(triplet);
      ++i;
    }
    else if (mid.jumps[j] < right.jumps[k]) {
      D[triplet] += mid.jumps[j] - prev_time;
      J[triplet] += 1.0;
      prev_time = mid.jumps[j];
      triplet = flip_mid_bit(triplet);
      ++j;
    }
    else {
      D[triplet] += right.jumps[k] - prev_time;
      J[triplet] += 1.0;
      prev_time = right.jumps[k];
      triplet = flip_right_bit(triplet);
      ++k;
    }
}


void
add_sufficient_statistics(const vector<Path> &paths,
                          vector<double> &J, vector<double> &D) {

  for (size_t i = 0; i < n_sites; ++i)
    paths[i].jumps.push_back(paths[i].tot_time);

  for (size_t i = 1; i < n_sites - 1; ++i)
    add_sufficient_statistics(paths[i-1], paths[i], paths[i+1], J, D);

  for (size_t i = 0; i < n_sites; ++i)
    paths[i].jumps.pop_back();
}

void
get_sufficient_statistics(const vector<vector<Path> > &the_paths,
                          vector<double> &J, vector<double> &D) {

  static const size_t n_triplets = 8;

  J.resize(n_triplets, 0.0);
  D.resize(n_triplets, 0.0);

  for (size_t i = 0; i < the_paths.size(); ++i)
    add_sufficient_statistics(the_paths[i], J, D);
}


string
format_two_by_two(const vector<vector<double> > &m) {
  std::ostringstream oss;
  oss << "["
      << std::setw(10) << std::right << m[0][0]
      << std::setw(10) << std::right << m[0][1]
      << "]\n["
      << std::setw(10) << std::right << m[1][0]
      << std::setw(10) << std::right << m[1][1]
      << "]";
  return oss.str();
}

// Steps:
// read in full paths
// extract initial seq, jumps
// sort jumps
// derive jump contexts

////////////////////////////////////////////////////////////////////////////////
// Estimate rates given branch lenghts and jumping times
////////////////////////////////////////////////////////////////////////////////

static double
log_likelihood(const vector<double> &J, const vector<double> &D,
               const vector<double> &rates) {

  static const size_t n_triplets = 8;

  assert(J.size() == n_triplets && D.size() == n_triplets &&
         rates.size() == n_triplets);

  double l = 0;
  for (size_t i = 0; i < n_triplets; ++i)
    l += J[i]*log(rates[i]) - D[i]*rates[i];

  return l;
}

/* compute gradients wrt log(rate[i]) */
void
get_gradient(const vector<double> &J, const vector<double> &D,
             const vector<double> &rates, vector<double> &gradient) {

  /* ordering of parameters:
     000 -> 010: 0 birth
     010 -> 000: 2 death
     001 -> 011: 1 expansion
     011 -> 001: 3 contraction
     100 -> 110: 4 expansion
     110 -> 100: 6 contraction
     101 -> 111: 5 merging
     111 -> 101: 7 splitting (not a free parameter)
   */
  static const size_t n_params = 8;
  gradient.resize(n_params, 0.0);

  const double factor_111 = (J[7] - D[7]*rates[7]);

  // [0] BIRTH: parameter 0 corresponds to 000->010
  gradient[0] = J[0] - D[0]*rates[0] + factor_111;

  // [2] DEATH: parameter 2 corresponds to 010->000
  gradient[2] = J[2] - D[2]*rates[2] - factor_111;

  // [1, 4] EXPANSION: parameter 1 corresponds to 001->011
  gradient[1] = J[1] + J[4] - (D[1] + D[4])*rates[1] - 2*factor_111;
  gradient[4] = gradient[1]; // expansion in other direction: 100->110

  // [3, 6] CONTRACTION: parameter 3 corresponds to 011->001
  gradient[3] = J[3] + J[6] - (D[3] + D[6])*rates[3] + 2*factor_111;
  gradient[6] = gradient[3]; // contraction in other direction: 110->100

  // [5] MERGING: parameter 2 corresponds to 101->111
  gradient[5] = J[5] - D[5]*rates[5] + factor_111;

  // [7] SPLITTING: parameter 7 corresponds to 111->101
  /* gradient[7] remains 0, since rates[7] is determined by other rates */
}

static void
candidate_rates(const double step,
                const vector<double> &rates,
                const vector<double> &gradient,
                vector<double> &updated_rates) {

  static const size_t n_rates = 8;

  updated_rates.resize(n_rates, 0);
  for (size_t i = 0; i < n_rates - 1; ++i)
    updated_rates[i] = exp(log(rates[i]) + gradient[i]*step);

  // final rate is in terms of other rates
  updated_rates[n_rates - 1] =
    exp(log(updated_rates[0])       // 000 -> 010 (once, numerator)
        + log(updated_rates[5])     // 101 -> 111 (once, numerator)
        + 2*log(updated_rates[3])   // 011 -> 001 (twice,numerator)
        - log(updated_rates[2])     // 010 -> 000 (once, denominator)
        - 2*log(updated_rates[1])); // 001 -> 011 (twice,denominator)
}

double
gradient_ascent(const double param_tol,
                const vector<double> &J,
                const vector<double> &D,
                const vector<double> &rates,
                vector<double> &updated_rates) {
  /* compute llk and gradient */
  double l = log_likelihood(J, D, rates);

  vector<double> g;
  get_gradient(J, D, rates, g);

  double norm = 0.0;
  for (size_t i = 0; i < g.size(); ++i)
    norm += fabs(g[i]);

  double step = 0.2/norm; // MAGIC!!!!
  double new_l = 0;
  do {
    step *= 0.5;
    candidate_rates(step, rates, g, updated_rates);
    new_l = log_likelihood(J, D, updated_rates);
  }
  while (new_l <= l && step > param_tol);

  return step;
}

void
estimate_rates(const double param_tol,
               const vector<vector<Path> > &the_paths,
               const vector<double> &rates,
               vector<double> &updated_rates) {

  vector<double> J;
  vector<double> D;
  get_sufficient_statistics(the_paths, J, D);

  vector<double> current_rates(rates);
  double current_llk = log_likelihood(J, D, current_rates);

  double diff = gradient_ascent(param_tol, J, D, current_rates, updated_rates);
  double updated_llk = log_likelihood(J, D, updated_rates);

  while (updated_llk - current_llk > param_tol && diff > param_tol) {
    current_llk = updated_llk;
    current_rates = updated_rates;
    diff = gradient_ascent(param_tol, J, D, current_rates, updated_rates);
    updated_llk = log_likelihood(J, D, updated_rates);
  }
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

static void
estimate_transition_probabilities(const vector<bool> &seq,
                                  vector<vector<double> > & trans_prob) {

  // count the different pair-wise transitions within the sequence
  double c00 = 0.0, c01 = 0.0, c10 = 0.0, c11 = 0.0;
  for (size_t i = 0; i < seq.size() - 1; ++i) {
    c00 += static_cast<size_t>((seq[i] == false) && (seq[i+1] == false));
    c01 += static_cast<size_t>((seq[i] == false) && (seq[i+1] == true));
    c10 += static_cast<size_t>((seq[i] == true) && (seq[i+1] == false));
    c11 += static_cast<size_t>((seq[i] == true) && (seq[i+1] == true));
  }
  if (c00*c11*c10*c01 == 0)
    cerr << "WARNING: Root sequence lack diversity" << endl;

  // calculate horizontal probability of no change
  const double p00 = c00/(c00 + c01);
  const double p11 = c11/(c10 + c11);

  // assign horizontal probabilities to matrix
  trans_prob = vector<vector<double> >(2, vector<double>(2, 0.0));
  trans_prob[0][0] = p00;
  trans_prob[0][1] = 1.0 - p00;
  trans_prob[1][0] = 1.0 - p11;
  trans_prob[1][1] = p11;
}

/* This function scales the branches of the tree so that the sum of
 * the branch lengths is 1.
 */
static void
scale_treesize(vector<double> &rates, vector<double> &branches) {

  // get the total sum of branch lengths in the tree
  const double treesize = accumulate(branches.begin(), branches.end(), 0.0);

  // scale the rates by the sum of the tree branches
  for (size_t i = 0; i < rates.size(); ++i)
    rates[i] *= treesize;

  // scale the tree branch lengths so they have unit sum
  for (size_t i = 0; i < branches.size(); ++i)
    branches[i] /= treesize;
}


int main(int argc, const char **argv) {

  try {

    string outfile;
    string param_file;
    bool VERBOSE = false;

    OptionParser opt_parse(strip_path(argv[0]), "estimate parameters from"
                           " complete data (site-specific paths)",
                           "<path-file>");
    opt_parse.add_opt("param", 'p', "params file", false, param_file);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
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
    if (leftover_args.size() < 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string path_file(leftover_args.back());
    ////////////////////////////////////////////////////////////////////////

    if (VERBOSE)
      cerr << "[READING PATHS: " << path_file << "]" << endl;
    vector<vector<Path> > paths; // along multiple branches
    vector<string> node_names;
    read_paths(path_file, node_names, paths);
    const size_t n_nodes = node_names.size();
    const size_t n_sites = paths.front().size();

    if (VERBOSE)
      cerr << "n_nodes=" << n_nodes << endl
           << "n_sites=" << n_sites << endl;

    vector<double> branches_from_paths;
    for (size_t i = 0; i < paths.size(); ++i)
      branches_from_paths.push_back(paths[i][0].tot_time);

    vector<bool> seq;
    get_initial_seq(paths[0], seq);

    vector<vector<double> > init_trans_prob;
    estimate_transition_probabilities(seq, init_trans_prob);

    if (VERBOSE)
      cerr << "horizontal transitions at root" << endl
           << format_two_by_two(init_trans_prob) << endl;

    vector<vector<Jump> > jumps;
    vector<Hold> holds;

    /* COLLECT ALL JUMPS */
    for (size_t node_id = 1; node_id < n_nodes; ++node_id) {

      vector<Jump> j;
      Hold h;
      get_jumps(paths[node_id], j, h);
      jumps.push_back(j);
      holds.push_back(h);

      if (VERBOSE) {
        cerr << "[PROCESSING BRANCH: " << node_names[node_id] << ", "
             << "SITES: " << n_sites << ", "
             << "JUMPS: " << jumps[node_id].size() << endl;
      }
    }

    /* IF MODEL PARAMETERS ARE PROVIDED, READ THEM */
    if (!param_file.empty()) {
      model_param p;
      p.read_param(param_file);
      vector<double> J;
      vector<double> D;
      get_sufficient_statistics(the_paths, J, D);

      vector<double> rates;
      p.get_rates(rates);

      vector<double> branches;
      p.t.get_branch_lengths(branches);
      branches.erase(branches.begin());

      if (VERBOSE) {
        cerr << "Provided starting rates:\t";
        for (size_t i = 0; i < rates.size(); ++i) cerr << rates[i] << "\t";
        cerr << endl;
        cerr << "Provided starting branches:\t";
        for (size_t i = 0; i < branches.size(); ++i) cerr << branches[i] << "\t";
        cerr << endl;
        cerr << "log-likelihood= "
             << log_likelihood(J, D, rates) << endl;
      }

      vector<double> gradient;
      get_gradient(J, D, rates, gradient);

      if (DEBUG) {
        cerr << "gradient :" << endl;
        for (size_t i = 0; i < gradient.size(); ++i)
          cerr << "[" << i << "]\t" << gradient[i] << endl;
      }

      double param_tol = 1e-10;
      vector<double> updated_rates;
      estimate_rates(param_tol, the_paths, rates, updated_rates);

      if (VERBOSE) {
        cerr << "optimized rates:\t" << endl;
        for (size_t i = 0; i < updated_rates.size(); ++i)
          cerr << updated_rates[i] << "\t";
        cerr << endl << "optimized likelihood:\t"
             << log_likelihood(J, D, updated_rates) << endl;
      }

      /* scale rates and branches to have unit branch length
         corresponding to one expected transition per site */
      const double the_rate_factor = rate_factor(updated_rates);
      vector<double> scaled_rates = updated_rates;
      vector<double> scaled_D = D;
      for (size_t i = 0; i < updated_rates.size(); ++i) {
        scaled_rates[i] /= the_rate_factor;
        scaled_D[i] *= the_rate_factor;
      }

      vector<double> scaled_branches(branches_from_paths);
      for (size_t b = 0; b < branches_from_paths.size(); ++b)
        scaled_branches[b] *= the_rate_factor;

      if (VERBOSE) {
        cerr << "scaled rates:\t" << endl;
        for (size_t i = 0; i < updated_rates.size(); ++i)
          cerr << scaled_rates[i] << "\t";
        cerr << endl;
        cerr << "scaled branches:\t" << endl;
        for (size_t i = 0; i < branches.size(); ++i)
          cerr << scaled_branches[i] << "\t";

        cerr << endl
             <<"new log-likelihood= "
             << log_likelihood(J, scaled_D, scaled_rates) << endl;
      }
    }
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
