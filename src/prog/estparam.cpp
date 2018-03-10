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

bool DEBUG = false;

// Steps:
// read in full paths
// extract initial seq, jumps
// sort jumps
// derive jump contexts

////////////////////////////////////////////////////////////////////////////////
// Estimate rates given branch lenghts and jumping times
////////////////////////////////////////////////////////////////////////////////

static double
log_likelihood(const vector<double> &J,
               const vector<double> &D,
               const vector<double> &rates) {

  static const size_t n_triplets = 8;

  assert(J.size() == n_triplets && D.size() == n_triplets &&
         rates.size() == n_triplets);

  double l = 0;
  for (size_t i = 0; i < n_triplets; ++i) {
    l += J[i]*log(rates[i]) - D[i]*rates[i];
  }
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
                vector<double> &new_rates) {

  static const size_t n_rates = 8;

  new_rates.resize(n_rates, 0);
  for (size_t i = 0; i < n_rates - 1; ++i)
    new_rates[i] = exp(log(rates[i]) + gradient[i]*step);

  // final rate is in terms of other rates
  new_rates[n_rates - 1] =
    exp(log(new_rates[0])       // 000 -> 010 (once, numerator)
        + log(new_rates[5])     // 101 -> 111 (once, numerator)
        + 2*log(new_rates[3])   // 011 -> 001 (twice,numerator)
        - log(new_rates[2])     // 010 -> 000 (once, denominator)
        - 2*log(new_rates[1])); // 001 -> 011 (twice,denominator)
}

double
gradient_ascent(const double param_tol,
                const vector<vector<Jump> > &jumps,
                const vector<double> &J,
                const vector<double> &D,
                const vector<double> &rates,
                vector<double> &new_rates) {
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
    candidate_rates(step, rates, g, new_rates);
    new_l = log_likelihood(J, D, new_rates);
  } while (new_l <= l && step > param_tol);

  return step;
}

void
estimate_rates(const double param_tol,
               const vector<vector<Jump> > &jumps,
               const vector<Hold> &holds,
               const vector<double> &rates,
               vector<double> &new_rates) {

  vector<double> J;
  vector<double> D;
  get_suff_stat(jumps, holds, J, D);

  vector<double> current_rates = rates;
  double current_llk = log_likelihood(J, D, current_rates);
  double diff = gradient_ascent(param_tol, jumps, J, D,
                                current_rates, new_rates);
  double new_llk = log_likelihood(J, D, new_rates);

  while (new_llk - current_llk > param_tol && diff > param_tol) {
    current_llk = new_llk;
    current_rates = new_rates;
    diff = gradient_ascent(param_tol, jumps, J, D, current_rates, new_rates);
    new_llk = log_likelihood(J, D, new_rates);

    if (DEBUG) {
      cerr << "###rates:\t";
      for (size_t i = 0; i < 8; ++i)
        cerr << new_rates[i] << "\t";
      cerr << new_llk << "\t" << new_llk - current_llk << endl;
    }
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

    OptionParser opt_parse(strip_path(argv[0]), "estimate parameter from paths",
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

    if (VERBOSE)
      cerr << "Reading paths " << path_file << endl;

    vector<vector<Path> > paths; // along multiple branches
    read_paths(path_file, paths);

    vector<double> branches_from_paths;
    for (size_t i = 0; i < paths.size(); ++i)
      branches_from_paths.push_back(paths[i][0].tot_time);

    vector<bool> seq;
    get_initial_seq(paths[0], seq);

    vector<vector<double> > trans_prob;
    estimate_transition_probabilities(seq, trans_prob);

    if (VERBOSE) {
      cerr << "Root transition probabilities are" << endl
           << "[";
      for (size_t i = 0; i < 2; ++i) {
        if (i > 0) cerr << " ";
        for (size_t j = 0; j < 2; ++j)
          cerr << trans_prob[i][j] << "\t";
        if (i > 0) cerr << "]";
        cerr << endl;
      }
    }

    vector<vector<Jump> > jumps;
    vector<Hold> holds;
    for (size_t b = 0; b < paths.size(); ++b) {
      vector<Jump> j;
      Hold h;
      get_jumps(paths[b], j, h);
      jumps.push_back(j);
      holds.push_back(h);

      if (VERBOSE) {
        cerr<< "Branch-" << b << endl;
        cerr << "Total sites: " << paths[b].size() << endl;
        cerr << "Total jumps: " << jumps[b].size() << endl;
      }

      //////////
      if (DEBUG) {
        cerr << "==========" << endl
             << "First 10 jumps:" << endl;
        for (size_t i = 0; i < 10; ++i)
          cerr << "position-" << jumps[b][i].pos
               << "\ttime-" << jumps[b][i].t
               << "\tcontext-" << jumps[b][i].context << endl;
        cerr << "==========" << endl;
        cerr << "==========" << endl
             << "Last 10 jumps:" << endl;
        for (size_t i = 10; i > 0; --i)
          cerr << "position-" << jumps[b][jumps[b].size()-i].pos
               << "\ttime-" << jumps[b][jumps[b].size()-i].t
               << "\tcontext-" << jumps[b][jumps[b].size()-i].context << endl;
        cerr << "==========" << endl;
      }
      //////////
    }

    if (!param_file.empty()) {
      model_param p;
      p.read_param(param_file);
      vector<double> J;
      vector<double> D;
      get_suff_stat(jumps, holds, J, D);

      //////////
      if (DEBUG) {
        cerr << "J:\t";
        size_t sum = 0;
        for (size_t i = 0; i < J.size(); ++i) {
          cerr << J[i] << "\t";
          sum += J[i];
        }
        cerr << "(" << sum << ")" << endl;

        for (size_t i = 0; i < 10; ++i) {
          cerr << "J" << i << "\tT(" << jumps[0][i].t << ")\t";
          for (size_t j = 0; j < 8; ++j)
            cerr << jumps[0][i].freq[j] << "\t";
          cerr << endl;
        }
        for (size_t i = 10; i > 0; --i) {
          cerr << "J" << jumps[0].size()-i
               << "\tT(" << jumps[0][jumps[0].size()-i].t << ")\t";
          for (size_t j = 0; j < 8; ++j)
            cerr << jumps[0][jumps[0].size()-i].freq[j] << "\t";
          cerr << endl;
        }
      }
      //////////

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
      if (DEBUG) {
        cerr << "tot_freq:" << endl;
        for (size_t i = 0; i < 8; ++i)
          cerr << "[" << i << "]\t" << J[i] << endl;
      }

      vector<double> gradient;
      get_gradient(J, D, rates, gradient);

      if (DEBUG) {
        cerr << "gradient :" << endl;
        for (size_t i = 0; i < gradient.size(); ++i)
          cerr << "[" << i << "]\t" << gradient[i] << endl;
      }

      double param_tol = 1e-10;
      vector<double> new_rates;
      estimate_rates(param_tol, jumps, holds, rates, new_rates);

      if (VERBOSE) {
        cerr << "optimized rates:\t" << endl;
        for (size_t i = 0; i < new_rates.size(); ++i)
          cerr << new_rates[i] << "\t";
        cerr << endl << "optimized likelihood:\t"
             << log_likelihood(J, D, new_rates) << endl;
      }

      /* scale rates and branches to have unit branch length
         corresponding to one expected transition per site */
      const double the_rate_factor = rate_factor(new_rates);
      vector<double> scaled_rates = new_rates;
      vector<double> scaled_D = D;
      for (size_t i = 0; i < new_rates.size(); ++i) {
        scaled_rates[i] /= the_rate_factor;
        scaled_D[i] *= the_rate_factor;
      }

      vector<double> scaled_branches(branches_from_paths);
      for (size_t b = 0; b < branches_from_paths.size(); ++b)
        scaled_branches[b] *= the_rate_factor;

      if (VERBOSE) {
        cerr << "scaled rates:\t" << endl;
        for (size_t i = 0; i < new_rates.size(); ++i)
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
