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
#include <cmath>   /* exp, sqrt, pow fabs*/
#include <numeric>  /* std::inner_product, accumulate*/
#include <functional>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "PhyloTreePreorder.hpp"
#include "Path.hpp"
#include "StateSeq.hpp"
#include "EpiEvoModel.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;


/* add the total time for each branch to the set of jump times at each
   site for each branch. This allows more convenient processing of
   jump sequences, but should be undone for cases when the final jump
   would imply a state change */
static void
add_total_time_to_jumps(vector<vector<Path> > &all_paths) {
  for (size_t i = 0; i < all_paths.size(); ++i)
    if (!all_paths[i].empty()) {
      // make sure all time durations are the same for each site
      assert(all_paths[i].front().tot_time ==
             all_paths[i].back().tot_time);
      const double tot_time = all_paths[i].front().tot_time;
      for (size_t j = 0; j < all_paths[i].size(); ++j)
        all_paths[i][j].jumps.push_back(tot_time);
    }
}

/* undo the changes made by "add_total_time_to_jumps" */
static void
remove_total_time_from_jumps(vector<vector<Path> > &all_paths) {
  for (size_t i = 0; i < all_paths.size(); ++i)
    if (!all_paths[i].empty()) {
      // make sure we are removing something added
      assert(all_paths[i].front().jumps.back() ==
             all_paths[i].back().jumps.back());
      for (size_t j = 0; j < all_paths[i].size(); ++j)
        all_paths[i][j].jumps.pop_back();
    }
}

void
add_sufficient_statistics(const Path &left, const Path &mid, const Path &right,
                          vector<double> &J, vector<double> &D) {

  // make sure the total time was put at the end of all jump sequences
  assert(left.jumps.back() == mid.jumps.back() &&
         left.jumps.back() == right.jumps.back());

  size_t triplet = triple2idx(left.init_state, mid.init_state, right.init_state);
  double prev_time = 0.0;
  size_t i = 0, j = 0, k = 0;
  while (i < left.jumps.size() && j < mid.jumps.size() && k < right.jumps.size())
    if (left.jumps[i] < std::min(mid.jumps[j], right.jumps[k])) { // LEFT
      D[triplet] += left.jumps[i] - prev_time;
      /* no need to update J[triplet] += 1.0; here */
      prev_time = left.jumps[i];
      triplet = flip_left_bit(triplet);
      ++i;
    }
    else if (mid.jumps[j] < right.jumps[k]) { // MID
      D[triplet] += mid.jumps[j] - prev_time;
      J[triplet] += 1.0;
      prev_time = mid.jumps[j];
      triplet = flip_mid_bit(triplet);
      ++j;
    }
    else { // RIGHT (also: default case, shouldn't udpate J)
      D[triplet] += right.jumps[k] - prev_time;
      /* no need to update J[triplet] += 1.0; here */
      prev_time = right.jumps[k];
      triplet = flip_right_bit(triplet);
      ++k;
    }
}

void
add_sufficient_statistics(const vector<Path> &paths,
                          vector<double> &J, vector<double> &D) {
  // iterate over sites with valid triples (i.e. not the first and last)
  const size_t n_sites = paths.size();
  assert(n_sites > 0); // make sure the root isn't processed
  for (size_t i = 1; i < n_sites - 1; ++i)
    add_sufficient_statistics(paths[i-1], paths[i], paths[i+1], J, D);
}


////////////////////////////////////////////////////////////////////////////////
// Estimate rates only
////////////////////////////////////////////////////////////////////////////////

void
get_sufficient_statistics(const vector<vector<Path> > &all_paths,
                          vector<double> &J, vector<double> &D) {

  static const size_t n_triplets = 8;

  J.resize(n_triplets, 0.0);
  D.resize(n_triplets, 0.0);

  // iterate over nodes, starting at 1 to avoid root
  for (size_t i = 1; i < all_paths.size(); ++i)
    add_sufficient_statistics(all_paths[i], J, D);
}


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


/* makes a single step of gradient ascent; identifies a new set of
   rates using the summary statistics in J and D */
double
gradient_ascent(const double param_tol,
                const vector<double> &J, const vector<double> &D,
                const vector<double> &rates, vector<double> &updated_rates) {

  /* compute llk and gradient */
  double llk = log_likelihood(J, D, rates);
  vector<double> gradient;
  get_gradient(J, D, rates, gradient);

  /* get the step size */
  double norm = 0.0;
  for (size_t i = 0; i < gradient.size(); ++i)
    norm += fabs(gradient[i]);
  double step_size = 0.2/norm; // MAGIC!!!!

  double updated_llk = 0;
  do {
    step_size *= 0.5;
    candidate_rates(step_size, rates, gradient, updated_rates);
    updated_llk = log_likelihood(J, D, updated_rates);
  }
  while (updated_llk <= llk && step_size > param_tol);

  return step_size;
}


double
estimate_rates(const double param_tol,
               const vector<double> &J, const vector<double> &D,
               const vector<double> &input_rates,
               vector<double> &rates) {

  rates = input_rates;
  double llk = log_likelihood(J, D, rates);

  bool improving = true;
  while (improving) {

    vector<double> updated_rates(rates);
    const double diff = gradient_ascent(param_tol, J, D, rates, updated_rates);
    const double updated_llk = log_likelihood(J, D, updated_rates);

    improving = (updated_llk - llk > param_tol && diff > param_tol);
    if (improving) {
      llk = updated_llk;
      rates = updated_rates;
    }
  }
  return llk;
}

/* This function scales the branches of the tree so that the sum of
 * the branch lengths is 1.
 */
static void
scale_treesize(vector<double> &rates, vector<double> &branches,
               vector<double> &D) {

  // get the total sum of branch lengths in the tree
  const double treesize = accumulate(branches.begin(), branches.end(), 0.0);

  // scale the rates by the sum of the tree branches
  for (size_t i = 0; i < rates.size(); ++i) {
    rates[i] *= treesize;
    D[i] /= treesize;
  }

  // scale the tree branch lengths so they have unit sum
  for (size_t i = 0; i < branches.size(); ++i)
    branches[i] /= treesize;
}


static void
scale_by_rate_factor(vector<double> &rates, vector<double> &branches,
                     vector<double> &D) {

  /* scale rates and branches to have unit branch length
     corresponding to one expected transition per site */
  const double the_rate_factor = rate_scaling_factor(rates);

  // scale the rates and summary stats for each triplet
  for (size_t triplet_id = 0; triplet_id < rates.size(); ++triplet_id) {
    rates[triplet_id] /= the_rate_factor;
    D[triplet_id] *= the_rate_factor;
  }

  // scale the length of the branch associated with each node
  for (size_t node_id = 0; node_id < branches.size(); ++node_id)
    branches[node_id] *= the_rate_factor;
}


////////////////////////////////////////////////////////////////////////////////
// Estimate rates and branches
////////////////////////////////////////////////////////////////////////////////

/* branches in likelihood and gradient ascent steps are scalers that should
   be multiplied to input branch lengths after the optimization is over */

void
get_sufficient_statistics(const vector<vector<Path> > &all_paths,
                          vector<vector<double> > &J,
                          vector<vector<double> > &D) {

  static const size_t n_triplets = 8;
  const size_t n_branches = all_paths.size();

  J.resize(n_branches, vector<double>(n_triplets, 0.0));
  D.resize(n_branches, vector<double>(n_triplets, 0.0));

  // iterate over nodes, starting at 1 to avoid root
  for (size_t i = 1; i < all_paths.size(); ++i)
    add_sufficient_statistics(all_paths[i], J[i], D[i]);
}


static double
log_likelihood(const vector<vector<double> > &J,
               const vector<vector<double> > &D,
               const vector<double> &rates,
               const vector<double> &branch_scalers) {

  static const size_t n_triplets = 8;
  const size_t n_branches = branch_scalers.size();

  assert(J.size() && J[0].size() == n_triplets &&
         D.size() && D[0].size() == n_triplets &&
         rates.size() == n_triplets);

  double l = 0;
  for (size_t b = 1; b < n_branches; ++b) {
    for (size_t i = 0; i < n_triplets; ++i) {
      l += (J[b][i]*log(rates[i]*branch_scalers[b]) -
            D[b][i]*rates[i]*branch_scalers[b]);
    }
  }

  return l;
}


/* compute gradients wrt log(rate[i]) and log(branch_scalers[b])*/
void
get_gradient(const vector<vector<double> > &J,
             const vector<vector<double> > &D,
             const vector<double> &rates,
             const vector<double> &branch_scalers,
             vector<double> &gradient) {

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

  for (size_t b = 1; b < branch_scalers.size(); ++b) {
    const double factor_111 = (J[b][7] - D[b][7]*rates[7]*branch_scalers[b]);

    // [0] BIRTH: parameter 0 corresponds to 000->010
    gradient[0] += J[b][0] - D[b][0]*rates[0]*branch_scalers[b] + factor_111;

    // [2] DEATH: parameter 2 corresponds to 010->000
    gradient[2] += J[b][2] - D[b][2]*rates[2]*branch_scalers[b] - factor_111;

    // [1, 4] EXPANSION: parameter 1 corresponds to 001->011
    gradient[1] += J[b][1] + J[b][4] - (D[b][1] + D[b][4])*rates[1]*branch_scalers[b] - 2*factor_111;
    gradient[4] = gradient[1]; // expansion in other direction: 100->110

    // [3, 6] CONTRACTION: parameter 3 corresponds to 011->001
    gradient[3] += J[b][3] + J[b][6] - (D[b][3] + D[b][6])*rates[3]*branch_scalers[b] + 2*factor_111;
    gradient[6] = gradient[3]; // contraction in other direction: 110->100

    // [5] MERGING: parameter 2 corresponds to 101->111
    gradient[5] += J[b][5] - D[b][5]*rates[5]*branch_scalers[b] + factor_111;

    // [7] SPLITTING: parameter 7 corresponds to 111->101
    /* gradient[7] remains 0, since rates[7] is determined by other rates */
  }
}


static void
candidate_rates(const double step,
                const vector<vector<double> > &J,
                const vector<vector<double> > &D,
                const vector<double> &rates,
                const vector<double> &branch_scalers,
                const vector<double> &gradient,
                vector<double> &updated_rates,
                vector<double> &updated_branch_scalers) {

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

  updated_branch_scalers = vector<double>(branch_scalers.size(), 0.0);

  // update branch_scalers with new rates
  for (size_t b = 1; b < branch_scalers.size(); ++b) {
    double numerator = 0.0;
    double denominator = 0.0;
    for (size_t i = 0; i < n_rates; ++i) {
      numerator += J[b][i];
      denominator += D[b][i]*updated_rates[i];
    }
    updated_branch_scalers[b] = numerator/denominator;
  }
}


/* makes a single step of gradient ascent; identifies a new set of
   rates and branch scalers using the summary statistics in J and D */
double
gradient_ascent(const double param_tol,
                const vector<vector<double> > &J,
                const vector<vector<double> > &D,
                const vector<double> &rates,
                const vector<double> &branch_scalers,
                vector<double> &updated_rates,
                vector<double> &updated_branch_scalers) {

  /* compute llk and gradient */
  double llk = log_likelihood(J, D, rates, branch_scalers);
  vector<double> gradient;
  get_gradient(J, D, rates, branch_scalers, gradient);

  /* get the step size */
  double norm = 0.0;
  for (size_t i = 0; i < gradient.size(); ++i)
    norm += fabs(gradient[i]);
  double step_size = 2.0/norm; // MAGIC!!!!

  double updated_llk = 0;
  do {
    step_size *= 0.5;
    candidate_rates(step_size, J, D, rates, branch_scalers, gradient,
                    updated_rates, updated_branch_scalers);
    updated_llk = log_likelihood(J, D, updated_rates, updated_branch_scalers);
  }
  while (updated_llk <= llk && step_size > param_tol);

  return step_size;
}


double
estimate_rates(const double param_tol,
               const vector<vector<double> > &J,
               const vector<vector<double> > &D,
               const vector<double> &input_rates,
               const vector<double> &input_branch_scalers,
               vector<double> &rates,
               vector<double> &branch_scalers) {

  rates = input_rates;
  branch_scalers = input_branch_scalers;
  double llk = log_likelihood(J, D, rates, branch_scalers);

  cerr << "llk=" << llk << endl;

  bool improving = true;
  while (improving) {

    vector<double> updated_rates(rates);
    vector<double> updated_branch_scalers(branch_scalers);
    const double diff = gradient_ascent(param_tol, J, D, rates, branch_scalers,
                                        updated_rates, updated_branch_scalers);
    const double updated_llk = log_likelihood(J, D, updated_rates, updated_branch_scalers);

    improving = (updated_llk - llk > param_tol && diff > param_tol);
    if (improving) {
      llk = updated_llk;
      rates = updated_rates;
      branch_scalers = updated_branch_scalers;
    }
  }

  cerr << "llk=" << llk << endl;

  return llk;
}

/* This function scales the branches of the tree so that the sum of
 * the branch lengths is 1.
 */
static void
scale_treesize(const vector<double> init_branches,
               vector<double> &rates, vector<double> &branches,
               vector<vector<double> > &D) {

  // get the total sum of branch lengths in the tree
  const double treesize = accumulate(branches.begin(), branches.end(), 0.0);

  // scale the rates by the sum of the tree branches
  for (size_t i = 0; i < rates.size(); ++i) {
    rates[i] *= treesize;
    for (size_t b = 0; b < D.size(); ++b)
      D[b][i] *= branches[b]/(init_branches[b]*treesize);
  }

  // scale the tree branch lengths so they have unit sum
  for (size_t i = 0; i < branches.size(); ++i)
    branches[i] /= treesize;
}

static void
scale_by_rate_factor(vector<double> &rates, vector<double> &branches,
                     vector<vector<double> > &D) {

  /* scale rates and branches to have unit branch length
     corresponding to one expected transition per site */
  const double the_rate_factor = rate_scaling_factor(rates);

  // scale the rates and summary stats for each triplet
  for (size_t triplet_id = 0; triplet_id < rates.size(); ++triplet_id) {
    rates[triplet_id] /= the_rate_factor;
    for (size_t b = 0; b < D.size(); ++b)
      D[b][triplet_id] *= the_rate_factor;
  }

  // scale the length of the branch associated with each node
  for (size_t node_id = 0; node_id < branches.size(); ++node_id)
    branches[node_id] *= the_rate_factor;
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



int main(int argc, const char **argv) {

  try {

    static const double param_tol = 1e-10;

    string outfile;
    bool VERBOSE = false;
    bool OPTBRANCH = false;

    ////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]), "estimate parameters from"
                           " complete data (site-specific paths)",
                           "<path-file> <init-param-file>");
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("branch", 'b', "optimize branch lengths as well",
                      false, OPTBRANCH);
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
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string path_file(leftover_args.front());
    const string init_model_file(leftover_args.back());
    ////////////////////////////////////////////////////////////////////////

    if (VERBOSE)
      cerr << "[READING PATHS: " << path_file << "]" << endl;
    vector<vector<Path> > all_paths; // along multiple branches
    vector<string> node_names;
    read_paths(path_file, node_names, all_paths);

    /* adding the total time to all jump sequences for convenience,
       but might need to remove later */
    add_total_time_to_jumps(all_paths);

    const size_t n_nodes = node_names.size();
    /* below: 1st element of all_paths empty at root; use last */
    const size_t n_sites = all_paths.back().size();

    if (VERBOSE)
      cerr << "n_nodes=" << n_nodes << endl
           << "n_sites=" << n_sites << endl;

    if (VERBOSE)
      cerr << "reading parameter file: " << init_model_file << endl;
    EpiEvoModel the_model;
    read_model(init_model_file, the_model);
    if (VERBOSE)
      cerr << the_model << endl;

    /* updated model branches with total_time in paths*/
    for (size_t i = 1; i < the_model.branches.size(); ++i) {
      the_model.branches[i] = all_paths[i].front().tot_time;
    }

    if (!OPTBRANCH) {

      if (VERBOSE)
        cerr << "[COMPUTING INITIAL SUFFICIENT STATISTICS]" << endl;

      // get initial values of sufficient statistics
      vector<double> J;
      vector<double> D;
      get_sufficient_statistics(all_paths, J, D);

      if (VERBOSE)
        cerr << "J:" << endl
             << triplet_info_to_string(J) << endl
             << "D:" << endl
             << triplet_info_to_string(D) << endl;

      // get branch lengths implicit in the path information
      vector<double> branches_from_paths(all_paths.size(), 0.0);
      for (size_t i = 0; i < all_paths.size(); ++i)
        if (!all_paths[i].empty())
          branches_from_paths[i] = all_paths[i][0].tot_time;

      if (VERBOSE)
        cerr << "initial log-likelihood: "
             << log_likelihood(J, D, the_model.triplet_rates) << endl;

      vector<double> gradient;
      get_gradient(J, D, the_model.triplet_rates, gradient);
      if (VERBOSE)
        cerr << "gradient:" << endl
             << triplet_info_to_string(gradient) << endl;

      vector<double> updated_rates;
      estimate_rates(param_tol, J, D, the_model.triplet_rates, updated_rates);
      if (VERBOSE)
        cerr << "updated rates:" << endl
             << triplet_info_to_string(updated_rates) << endl;

      if (VERBOSE)
        cerr << "updated log-likelihood: "
             << log_likelihood(J, D, updated_rates) << endl;

      /* scale rates and branches to have unit branch length
         corresponding to one expected transition per site */
      vector<double> updated_branches(the_model.branches);
      // scale_by_rate_factor(updated_rates, updated_branches, D);
      scale_treesize(updated_rates, updated_branches, D);

      if (VERBOSE) {
        cerr << "scaled branches:\t" << endl;
        for (size_t i = 0; i < node_names.size(); ++i)
          cerr << node_names[i] << '\t' << updated_branches[i] << endl;

        cerr << "scaled rates:\t" << endl
             << triplet_info_to_string(updated_rates) << endl;

        cerr << "scaled D:" << endl
             << triplet_info_to_string(D) << endl;

        if (VERBOSE)
          cerr << "scaled log-likelihood: "
               << log_likelihood(J, D, updated_rates) << endl;
      }

    } else {

      if (VERBOSE)
        cerr << "[COMPUTING INITIAL SUFFICIENT STATISTICS]" << endl;

      vector<vector<double> > J;
      vector<vector<double> > D;
      get_sufficient_statistics(all_paths, J, D);

      if (VERBOSE)
        for (size_t b = 1; b < all_paths.size(); ++b)
          cerr << "J["<< b << "]:" << endl
               << triplet_info_to_string(J[b]) << endl
               << "D["<< b << "]:" << endl
               << triplet_info_to_string(D[b]) << endl;

      vector<double> gradient;
      get_gradient(J, D, the_model.triplet_rates,
                   the_model.branches, gradient);
      if (VERBOSE)
        cerr << "gradient:" << endl
             << triplet_info_to_string(gradient) << endl;

      vector<double> updated_rates;
      const vector<double> initial_branches(the_model.branches);
      vector<double> initial_branch_scalers(the_model.branches.size(), 1.0);
      vector<double> updated_branch_scalers;

      estimate_rates(param_tol, J, D, the_model.triplet_rates,
                     the_model.branches, updated_rates, updated_branch_scalers);

      for (size_t i = 0; i < node_names.size(); ++i)
        the_model.branches[i] *= updated_branch_scalers[i];

      if (VERBOSE) {
        cerr << "updated rates:" << endl
             << triplet_info_to_string(updated_rates) << endl;
        cerr << "updated_branch_scalers:\t" << endl;
        for (size_t i = 0; i < node_names.size(); ++i)
          cerr << node_names[i] << '\t' << updated_branch_scalers[i] << endl;

        cerr << "updated_branches:\t" << endl;
        for (size_t i = 0; i < node_names.size(); ++i)
          cerr << node_names[i] << '\t' << the_model.branches[i] << endl;
      }


      if (VERBOSE)
        cerr << "updated log-likelihood: "
             << log_likelihood(J, D, updated_rates, updated_branch_scalers) << endl;

      /* scale rates and branches to have unit branch length
         corresponding to one expected transition per site */
      vector<double> updated_branches(the_model.branches);
      // scale_by_rate_factor(updated_rates, updated_branches, D);
      scale_treesize(initial_branches, updated_rates, updated_branches, D);

      if (VERBOSE) {
        cerr << "scaled branches:\t" << endl;
        for (size_t i = 0; i < node_names.size(); ++i)
          cerr << node_names[i] << '\t' << updated_branches[i] << endl;
        cerr << "scaled rates:\t" << endl
             << triplet_info_to_string(updated_rates) << endl;

        // cerr << "scaled D:" << endl;
        // for (size_t b = 1; b < D.size(); ++b )
        //   cerr << "branch " << node_names[b] << endl
        //        << triplet_info_to_string(D[b]) << endl;
        if (VERBOSE)
          cerr << "scaled log-likelihood: "
               << log_likelihood(J, D, updated_rates, initial_branch_scalers) << endl;
      }

    }

  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
