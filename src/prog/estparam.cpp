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

double llk(const vector<double> &J,
           const vector<double> &D,
           const vector<double> &rates) {
  double l = 0;
  for (size_t i = 0; i < 8; ++i) {
    l += J[i]*log(rates[i]) - D[i]*rates[i];
  }
  return l;
}

void get_gradient(const vector<double> &J,
                  const vector<double> &D,
                  const vector<double> &rates,
                  vector<double> &gradient) {
  /* gradients w.r.t log(rate[i])*/
  gradient.resize(8, 0);
  gradient[0] = J[0] - D[0]*rates[0] + J[7] - D[7]*rates[7];
  gradient[2] = J[2] - D[2]*rates[2] - J[7] + D[7]*rates[7];
  gradient[5] = J[5] - D[5]*rates[5] + J[7] - D[7]*rates[7];
  gradient[1] = J[1] + J[4] - (D[1]+D[4])*rates[1] - 2*J[7] + 2*D[7]*rates[7];
  gradient[3] = J[3] + J[6] - (D[3]+D[6])*rates[3] + 2*J[7] - 2*D[7]*rates[7];
  gradient[4] = gradient[1];
  gradient[6] = gradient[3];
  /* gradient[7] remains 0, since rates[7] is determined by other rates. */
}

void candidate_rates(const double param_tol, const double step,
                     const vector<double> &rates,
                     const vector<double> &gradient,
                     vector<double> &new_rates) {
  new_rates.resize(8, 0);
  for (size_t i = 0; i < 7; ++i)
    new_rates[i] = exp(log(rates[i]) + gradient[i]*step);

  new_rates[7] = exp(log(new_rates[0]) + log(new_rates[5]) +  2*log(new_rates[3])
                     - log(new_rates[2]) - 2*log(new_rates[1]));
}

double
gradient_ascent(const double param_tol,
                const vector<vector<Jump> > &jumps,
                const vector<double> &J,
                const vector<double> &D,
                const vector<double> &rates,
                vector<double> &new_rates) {
  /* compute llk and gradient */
  double l = llk(J, D, rates);

  vector<double> g;
  get_gradient(J, D, rates, g);
  double norm = 0.0;
  for (size_t i = 0; i < g.size(); ++i)
    norm += fabs(g[i]);

  double step = 0.2/norm;
  double new_l = 0;
  do {
    step *= 0.5;
    candidate_rates(param_tol, step, rates, g, new_rates);
    new_l = llk(J, D, new_rates);
  } while (new_l <= l && step > param_tol);
  return step;
}

void est_rates(const double param_tol,
               const vector<vector<Jump> > &jumps,
               const vector<Hold> &holds,
               const vector<double> &rates,
               vector<double> &new_rates) {
  vector<double> J;
  vector<double> D;
  get_suff_stat(jumps, holds, J, D);

  vector<double> current_rates = rates;
  double current_llk = llk(J, D, current_rates);
  double diff = gradient_ascent(param_tol, jumps, J, D,
                                current_rates, new_rates);
  double new_llk = llk(J, D, new_rates);

  while (new_llk - current_llk > param_tol && diff > param_tol) {
    current_llk = new_llk;
    current_rates = new_rates;
    diff = gradient_ascent(param_tol, jumps, J, D,
                           current_rates, new_rates);
    new_llk = llk(J, D, new_rates);

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

void est_trans_prob(const vector<bool> &seq,
                    vector<vector<double> > & trans_prob) {
  double c00(0), c01(0), c10(0), c11(0), c0(0), c1(1);
  for (size_t i = 0; i < seq.size() -1; ++i) {
    c00 += (size_t)((!seq[i]) & (!seq[i+1]));
    c01 += (size_t)((!seq[i]) & seq[i+1]);
    c10 += (size_t)(seq[i] & (!seq[i+1]));
    c11 += (size_t)(seq[i] & seq[i+1]);
    c0 += (size_t)(!seq[i]);
    c1 += (size_t)(seq[i]);
  }

  if (c00*c11*c10*c01 == 0) {
    cerr << "WARNING: Root sequence lack diversity" << endl;
  }

  const double p00 = c00/c0;
  const double p11 = c11/c1;
  trans_prob = vector<vector<double> >(2, vector<double>(2, 0.0));
  trans_prob[0][0] = p00;
  trans_prob[0][1] = 1 - p00;
  trans_prob[1][0] = 1 - p11;
  trans_prob[1][1] = p11;
}

void scale_treesize(vector<double> &rates, vector<double> &branches) {
  double treesize = std::accumulate(branches.begin(), branches.end(), 0.0);
  for (size_t i = 0; i < rates.size(); ++i)
    rates[i] *= treesize;
  for (size_t b = 0; b < branches.size(); ++b)
    branches[b] /= treesize;
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
    for (size_t i = 0; i < paths.size(); ++i) {
      branches_from_paths.push_back(paths[i][0].tot_time);
    }

    vector<bool> seq;
    get_inital_seq(paths[0], seq);

    vector<vector<double> > trans_prob;
    est_trans_prob(seq, trans_prob);

    if (VERBOSE) {
      cerr << "Root transition probabilities are" << endl
           << "[";
      for (size_t i = 0; i < 2; ++i) {
        if (i > 0) cerr << " ";
        for (size_t j = 0; j < 2; ++j) {
          cerr << trans_prob[i][j] << "\t";
        }
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
             << llk(J, D, rates) << endl;
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
      est_rates(param_tol, jumps, holds, rates, new_rates);

      if (VERBOSE) {
        cerr << "optimized rates:\t" << endl;
        for (size_t i = 0; i < new_rates.size(); ++i)
          cerr << new_rates[i] << "\t";
        cerr << endl << "optimized likelihood:\t"
             << llk(J, D, new_rates) << endl;
      }

      /* scale rates and branches to have unit branch length corresponding
         to one expected transition per site */
      double rf = rate_factor(new_rates);
      vector<double> scaled_rates = new_rates;
      vector<double> scaled_D = D;
      for (size_t i = 0; i < new_rates.size(); ++i) {
        scaled_rates[i] /= rf;
        scaled_D[i] *= rf;
      }

      vector<double> scaled_branches = branches_from_paths;
      for (size_t b = 0; b < branches_from_paths.size(); ++b)
        scaled_branches[b] *= rf;

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
             << llk(J, scaled_D, scaled_rates) << endl;
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
