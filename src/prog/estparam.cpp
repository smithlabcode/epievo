
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <math.h>   /* exp, sqrt, pow fabs*/
#include <numeric>  /* std::inner_product */
#include <algorithm>    // std::sort

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "PhyloTreePreorder.hpp"
#include "path.hpp"  /* related to Path */
#include "TwoStateSeq.hpp"
#include "param.hpp" /* model_param */
#include "TripletPattern.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;

bool DEBUG = true;

struct Jump {
  double t;
  size_t pos;
  size_t context;
  vector<size_t> freq;
  Jump(double time, size_t position): t(time), pos(position) {};
};

bool jump_compare(Jump lhs, Jump rhs) { return lhs.t < rhs.t;}

// Steps:
// read in full paths
// extract initial seq, jumps
// sort jumps
// derive jump contexts

void get_jumps(const vector<Path> &paths, vector<Jump> &jumps) {
  for (size_t i = 0; i < paths.size(); ++i) {
    for (size_t j = 0; j < paths[i].jumps.size(); ++j) {
      jumps.push_back(Jump(paths[i].jumps[j], i));
    }
  }
  std::sort(jumps.begin(), jumps.end(), jump_compare);


  // mutate from start and compute context and freqs at jumps
  vector<bool> seq;
  get_inital_seq(paths, seq);

  PatSeq patseq(seq);
  for (size_t i = 0; i < jumps.size(); ++i) {
    jumps[i].context = patseq.get_context(jumps[i].pos);
    patseq.get_all_context_freq(jumps[i].freq);
    patseq.mutate(jumps[i].pos);
  }

  // if (DEBUG) {
  //   vector<size_t> tf(8,0);
  //   for (size_t i = 2; i <seq.size(); ++i)
  //     ++tf[(size_t)(seq[i-2])*4 + (size_t)(seq[i-1])*2 + (size_t)(seq[i])];

  //   cerr << "Seq-freq:";
  //   for (size_t i = 0; i <tf.size(); ++i)
  //     cerr << tf[i] << "\t";
  //   cerr << endl;
  // }

  // TwoStateSeq tss(seq);
  // for (size_t i = 0; i < jumps.size(); ++i) {
  //   tss.count_s3(jumps[i].freq); // get context frequencies before the jump
  //   tss.mutate(jumps[i].pos, jumps[i].context); //get mutation context

  //   if (i==0 && DEBUG) {
  //     cerr << "Cnt-freq:";
  //     for (size_t j = 0; j <8; ++j)
  //       cerr << jumps[i].freq[j] << "\t";
  //     cerr << endl;
  //   }
  //
  // }

}

void get_suff_stat(const vector<vector<Jump> > &jumps,
                   vector<double> &tot_freq,
                   vector<double> &weights) {
  tot_freq.resize(8, 0);
  weights.resize(8, 0);
  double prev_jump_time = 0;
  for (size_t branch = 0; branch < jumps.size(); ++branch) {
    for (size_t i = 0; i < jumps[branch].size(); ++i) {
      ++tot_freq[jumps[branch][i].context];
      weights[jumps[branch][i].context] +=
        jumps[branch][i].freq[jumps[branch][i].context] *
        (jumps[branch][i].t - prev_jump_time);
      prev_jump_time = jumps[branch][i].t;
    }
  }
}


void get_suff_stat(const vector<Jump> &jumps,
                   vector<double> &tot_freq,
                   vector<double> &weights) {
  tot_freq.resize(8, 0);
  weights.resize(8, 0);
  double prev_jump_time = 0;
  for (size_t i = 0; i < jumps.size(); ++i) {
    ++tot_freq[jumps[i].context];
    weights[jumps[i].context] +=
      jumps[i].freq[jumps[i].context] * (jumps[i].t - prev_jump_time);
    prev_jump_time = jumps[i].t;
  }
}


double llk(const vector<vector<Jump> > &jumps,
           const vector<double> &tot_freq,
           const vector<double> &weights,
           const vector<double> &rates) {
  double l = 0;
  for (size_t branch = 0; branch < jumps.size(); ++branch) {
    for (size_t i = 0; i < 8; ++i) {
      l += tot_freq[i]*(log(rates[i]) -
                        log(jumps[branch][i].freq[jumps[branch][i].context])) -
        weights[i]*rates[i];
    }
  }

  return l;
}


double llk(const vector<Jump> &jumps,
           const vector<double> &tot_freq,
           const vector<double> &weights,
           const vector<double> &rates) {
  double l = 0;
  for (size_t i = 0; i < 8; ++i) {
    l += tot_freq[i]*(log(rates[i]) - log(jumps[i].freq[jumps[i].context])) -
      weights[i]*rates[i];
  }

  return l;
}


void get_gradient(const vector<double> &tot_freq,
                  const vector<double> &weights,
                  const vector<double> &rates,
                  vector<double> &gradient) {
  gradient.resize(8, 0);

  gradient[0] = (tot_freq[0] + tot_freq[7])/rates[0] - weights[0] -
    weights[7]*rates[7]/rates[0];

  gradient[2] = (tot_freq[2] - tot_freq[7])/rates[2] - weights[2] +
    weights[7]*rates[7]/rates[2];

  gradient[5] = (tot_freq[5] + tot_freq[7])/rates[5] - weights[5] -
    weights[7]*rates[7]/rates[5];

  gradient[1] = (tot_freq[1] + tot_freq[4] - 2*tot_freq[7])/rates[1] -
    (weights[1] + weights[4]) + 2*weights[7]*rates[7]/rates[1];

  gradient[3] = (tot_freq[3] + tot_freq[6] + 2*tot_freq[7])/rates[3] -
    (weights[3] + weights[6]) - 2*weights[7]*rates[7]/rates[3];

  gradient[4] = gradient[1];
  gradient[6] = gradient[3];
  // gradient[7] remains 0, since rates[7] is determined by other rates.
}


void candidate_rates(const double param_tol,
                     const double init_step,
                     const vector<double> &rates,
                     const vector<double> &gradient,
                     vector<double> &new_rates,
                     double &abs_diff) {
  new_rates.resize(rates.size(), 0);
  double norm = 0.0;
  for (size_t i = 0; i < gradient.size(); ++i) {
    norm += fabs(gradient[i]);
  }
  double step = init_step/norm;
  bool found = false;
  while (!found && step < param_tol) {
    step *= 0.5;
    for (size_t i = 0; i < rates.size() - 1; ++i) {
      new_rates[i] = rates[i] + gradient[i]*step;
    }
    new_rates.back() = new_rates[0]*new_rates[5]*
      pow(new_rates[3]/ new_rates[1], 2)/new_rates[2];
    found = true;
    for (size_t i = 0; i < rates.size(); ++i) {
      if (new_rates[i] <= 0) found = false;
    }
  }
  abs_diff = step*norm;
}


double gradient_ascent(const double param_tol,
                       const vector<vector<Jump> > &jumps,
                       const vector<double> &tot_freq,
                       const vector<double> &weights,
                       const vector<double> &rates,
                       vector<double> &new_rates) {
  // compute llk and gradient
  vector<double> g;
  double l = llk(jumps, tot_freq, weights, rates);
  get_gradient(tot_freq, weights, rates, g);
  double step = 0.2;
  double diff = 1;
  double new_l = 0;
  do {
    step *= 0.5;
    candidate_rates(param_tol, step, rates, g, new_rates, diff);
    new_l = llk(jumps, tot_freq, weights, new_rates);
  } while (new_l <= l && diff > param_tol);
  return diff;
}


double gradient_ascent(const double param_tol,
                       const vector<Jump> &jumps,
                       const vector<double> &tot_freq,
                       const vector<double> &weights,
                       const vector<double> &rates,
                       vector<double> &new_rates) {
  // compute llk and gradient
  vector<double> g;
  double l = llk(jumps, tot_freq, weights, rates);
  get_gradient(tot_freq, weights, rates, g);
  double step = 0.2;
  double diff = 1;
  double new_l = 0;
  do {
    step *= 0.5;
    candidate_rates(param_tol, step, rates, g, new_rates, diff);
    new_l = llk(jumps, tot_freq, weights, new_rates);
  } while (new_l <= l && diff > param_tol);
  return diff;
}


void est_param(const double param_tol,
               const vector<vector<Jump> > &jumps,
               const vector<double> &rates,
               vector<double> &new_rates) {
  vector<double> tot_freq;
  vector<double> weights;
  get_suff_stat(jumps, tot_freq, weights);

  vector<double> current_rates = rates;
  double current_llk = llk(jumps, tot_freq, weights, current_rates);
  double diff = gradient_ascent(param_tol, jumps, tot_freq,
                                weights, current_rates, new_rates);
  double new_llk = llk(jumps, tot_freq, weights, new_rates);



  while (new_llk - current_llk > param_tol && diff > param_tol) {
    current_llk = new_llk;
    current_rates = new_rates;
    diff = gradient_ascent(param_tol, jumps, tot_freq,
                           weights, current_rates, new_rates);
    new_llk = llk(jumps, tot_freq, weights, new_rates);
    cerr << "***rates:\t";
    for (size_t i = 0; i < new_rates.size(); ++i)
      cerr << new_rates[i] << "\t";
    cerr << new_llk << endl;
  }
}

void est_param(const double param_tol,
               const vector<Jump> &jumps,
               const vector<double> &rates,
               vector<double> &new_rates) {
  vector<double> tot_freq;
  vector<double> weights;
  get_suff_stat(jumps, tot_freq, weights);

  vector<double> current_rates = rates;
  double current_llk = llk(jumps, tot_freq, weights, current_rates);
  double diff = gradient_ascent(param_tol, jumps, tot_freq,
                                weights, current_rates, new_rates);
  double new_llk = llk(jumps, tot_freq, weights, new_rates);
  while (new_llk - current_llk > param_tol && diff > param_tol) {
    current_llk = new_llk;
    current_rates = new_rates;
    diff = gradient_ascent(param_tol, jumps, tot_freq,
                           weights, current_rates, new_rates);
    new_llk = llk(jumps, tot_freq, weights, new_rates);
  }
}

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
  const double p00 = c00/c0;
  const double p11 = c11/c1;
  trans_prob = vector<vector<double> >(2, vector<double>(2, 0.0));
  trans_prob[0][0] = p00;
  trans_prob[0][1] = 1 - p00;
  trans_prob[1][0] = 1 - p11;
  trans_prob[1][1] = p11;
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

    // const size_t n_sites = paths[0].size();
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
    for (size_t branch = 0; branch < paths.size(); ++branch) {
      vector<Jump> j;
      get_jumps(paths[branch], j);
      jumps.push_back(j);

      if (VERBOSE) {
      cerr << "Total sites: " << paths[branch].size() << endl;
      cerr << "Total jumps: " << jumps[branch].size() << endl;
      }
      if (DEBUG) {
        cerr << "==========" << endl
             << "First 10 jumps:" << endl;
        for (size_t i = 0; i < 10; ++i)
          cerr << "position-" << jumps[branch][i].pos
               << "\ttime-" << jumps[branch][i].t
               << "\tcontext-" << jumps[branch][i].context << endl;
        cerr << "==========" << endl;
      }
    }

    if (!param_file.empty()) {
      model_param p;
      read_param(param_file, p);
      vector<double> tot_freq;
      vector<double> weights;
      get_suff_stat(jumps, tot_freq, weights);

      if (DEBUG) {
        for (size_t i = 0; i < 10; ++i) {
          cerr << "J" << i << "\tT(" << jumps[0][i].t << ")\t";
          for (size_t j = 0; j < 8; ++j)
            cerr << jumps[0][i].freq[j] << "\t";
          cerr << endl;
        }
        for (size_t i = 0; i < 10; ++i) {
          cerr << "J" << jumps[0].size()-11+i << "\tT(" << jumps[0][i].t << ")\t";
          for (size_t j = 0; j < 8; ++j)
            cerr << jumps[0][jumps[0].size()-11+i].freq[j] << "\t";
          cerr << endl;
        }
      }

      vector<double> rates;
      get_rates(p, rates);

      if (DEBUG) {
        cerr << "rates:\t";
        for (size_t i = 0; i < rates.size(); ++i) cerr << rates[i] << "\t";
        cerr << endl;
        cerr << "log-likelihood= " << llk(jumps, tot_freq, weights, rates) << endl;

        cerr << "tot_freq:" << endl;
        for (size_t i =0; i < 8; ++i)
          cerr << "[" << i << "]\t" << tot_freq[i] << endl;

        cerr << "weights:" << endl;
        for (size_t i =0; i < 8; ++i)
          cerr << "[" << i << "]\t" << weights[i] << endl;

        vector<double> gradient;
        get_gradient(tot_freq, weights, rates, gradient);
        cerr << "gradient :" << endl;
        for (size_t i =0; i < gradient.size(); ++i)
          cerr << "[" << i << "]\t" << gradient[i] << endl;

        double param_tol = 1e-5;
        vector<double> new_rates;
        est_param(param_tol, jumps, rates, new_rates);
        cerr << "new rates:\t" << endl;
        for (size_t i = 0; i < new_rates.size(); ++i)
          cerr << new_rates[i] << "\t";
        cerr << endl
             <<"new log-likelihood= "
             << llk(jumps, tot_freq, weights, new_rates) << endl;
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
