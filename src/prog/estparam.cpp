
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <math.h>   /* exp, sqrt, pow fabs*/
#include <numeric>  /* std::inner_product */
#include <algorithm>    // std::sort
#include <functional>

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

bool DEBUG = false;

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
  /*only set jumping times and positions */
  for (size_t i = 0; i < paths.size(); ++i) {
    for (size_t j = 0; j < paths[i].jumps.size(); ++j) {
      jumps.push_back(Jump(paths[i].jumps[j], i));
    }
  }
  std::sort(jumps.begin(), jumps.end(), jump_compare);

  /* mutate from start and compute context and freqs at jumps */
  vector<bool> seq;
  get_inital_seq(paths, seq);

  PatSeq patseq(seq);
  for (size_t i = 0; i < jumps.size(); ++i) {
    jumps[i].context = patseq.get_context(jumps[i].pos);
    patseq.get_all_context_freq(jumps[i].freq);
    patseq.mutate(jumps[i].pos);
  }
}

void get_suff_stat(const vector<vector<Jump> > &jumps,
                   vector<double> &J,
                   vector<double> &D) {
  /* J_{ijk} = Total number of jumps in context ijk*/
  /* D_{ijk} = Sum over all jumps:
     (holding time before such a jump)*(frequency of jump context) */
  J.resize(8, 0);
  D.resize(8, 0);
  for (size_t branch = 0; branch < jumps.size(); ++branch) {
    double prev_jump_time = 0;
    const size_t n_jumps = jumps[branch].size();
    for (size_t i = 0; i < n_jumps; ++i) {
      size_t context = jumps[branch][i].context;
      vector<size_t> freq = jumps[branch][i].freq;
      double jump_time = jumps[branch][i].t ;
      if (jump_time - prev_jump_time > 0) {
        J[context] += 1.0;
        for (size_t ct = 0; ct < 8; ++ct) {
          D[ct] +=  freq[ct]*jump_time - freq[ct]*prev_jump_time;
        }
      }
      prev_jump_time = jump_time;
    }
  }
}


void get_suff_stat(const vector<Jump> &jumps,
                   vector<double> &J,
                   vector<double> &D) {
  J.resize(8, 0);
  D.resize(8, 0);
  double prev_jump_time = 0;
  for (size_t i = 0; i < jumps.size(); ++i) {
    size_t context = jumps[i].context;
    vector<size_t> freq = jumps[i].freq;
    double jump_time = jumps[i].t ;
    J[context] += 1 ;
    for (size_t ct = 0; ct < 8; ++ct) {
      D[ct] += freq[ct]*jump_time - freq[ct]*prev_jump_time;
    }
    prev_jump_time = jump_time;
  }
}

////////// likelihood w.r.t logrates //////////
double llk(const vector<double> &J,
           const vector<double> &D,
           const vector<double> &rates) {
  double l = 0;
  for (size_t i = 0; i < 8; ++i) {
    l += J[i]*log(rates[i]) - D[i]*rates[i];
  }
  return l;
}

void get_gradient_logrates(const vector<double> &J,
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
  // gradient[7] remains 0, since rates[7] is determined by other rates.
}

void candidate_logrates(const double param_tol,
                        const double step,
                        const vector<double> &rates,
                        const vector<double> &gradient,
                        vector<double> &new_rates) {
  new_rates.resize(8, 0);
  for (size_t i = 0; i < 7; ++i)
    new_rates[i] = exp(log(rates[i]) + gradient[i]*step);

  new_rates[7] = exp(log(new_rates[0]) + log(new_rates[5]) +  2*log(new_rates[3])
                     - log(new_rates[2]) - 2*log(new_rates[1]));
}

double gradient_ascent_logrates(const double param_tol,
                                const vector<vector<Jump> > &jumps,
                                const vector<double> &J,
                                const vector<double> &D,
                                const vector<double> &rates,
                                vector<double> &new_rates) {
  /* compute llk and gradient */
  double l = llk(J, D, rates);

  vector<double> g;
  get_gradient_logrates(J, D, rates, g);
  double norm = 0.0;
  for (size_t i = 0; i < g.size(); ++i)
    norm += fabs(g[i]);

  double step = 0.2/norm;
  double diff = 1;
  double new_l = 0;
  do {
    step *= 0.5;
    candidate_logrates(param_tol, step, rates, g, new_rates);
    new_l = llk(J, D, new_rates);
  } while (new_l <= l && diff > param_tol);
  return diff;
}


double gradient_ascent_logrates(const double param_tol,
                                const vector<Jump> &jumps,
                                const vector<double> &J,
                                const vector<double> &D,
                                const vector<double> &rates,
                                vector<double> &new_rates) {
  /* compute llk and gradient */
  vector<double> g;
  double l = llk(J, D, rates);
  get_gradient_logrates(J, D, rates, g);
  double step = 0.2;
  double new_l = 0;
  do {
    step *= 0.5;
    candidate_logrates(param_tol, step, rates, g, new_rates);
    new_l = llk(J, D, new_rates);
  } while (new_l <= l && step > param_tol);
  return step;
}

void est_param_logrates(const double param_tol,
                        const vector<vector<Jump> > &jumps,
                        const vector<double> &rates,
                        vector<double> &new_rates) {

  vector<double> J;
  vector<double> D;
  get_suff_stat(jumps, J, D);

  vector<double> current_rates = rates;
  double current_llk = llk(J, D, current_rates);
  double diff = gradient_ascent_logrates(param_tol, jumps, J, D,
                                         current_rates, new_rates);
  double new_llk = llk(J, D, new_rates);

  while (new_llk - current_llk > param_tol && diff > param_tol) {
    current_llk = new_llk;
    current_rates = new_rates;
    diff = gradient_ascent_logrates(param_tol, jumps, J, D,
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


// ////////// likelihood w.r.t rates //////////

// double llk(const vector<double> &J,
//            const vector<double> &D,
//            const vector<double> &rates) {
//   double l = -std::inner_product(D.begin(), D.end(), rates.begin(), 0);
//   for (size_t i = 0; i < 8; ++i) {
//     l += J[i]*log(rates[i]);
//   }
//   return l;
// }

// void get_gradient(const vector<double> &J,
//                   const vector<double> &D,
//                   const vector<double> &rates,
//                   vector<double> &gradient) {
//   gradient.resize(8, 0);
//   gradient[0] = (J[0] + J[7])/rates[0] - D[0] - D[7]*rates[7]/rates[0];
//   gradient[2] = (J[2] - J[7])/rates[2] - D[2] + D[7]*rates[7]/rates[2];
//   gradient[5] = (J[5] + J[7])/rates[5] - D[5] - D[7]*rates[7]/rates[5];
//   gradient[1] = (J[1] + J[4] - 2*J[7])/rates[1] - (D[1] + D[4]) +
//     2*D[7]*rates[7]/rates[1];
//   gradient[3] = (J[3] + J[6] + 2*J[7])/rates[3] - (D[3] + D[6]) -
//     2*D[7]*rates[7]/rates[3];
//   gradient[4] = gradient[1];
//   gradient[6] = gradient[3];
//   // gradient[7] remains 0, since rates[7] is determined by other rates.
// }


// void candidate_rates(const double param_tol,
//                      const double init_step,
//                      const vector<double> &rates,
//                      const vector<double> &gradient,
//                      vector<double> &new_rates,
//                      double &abs_diff) {
//   new_rates.resize(rates.size(), 0);
//   double norm = 0.0;
//   for (size_t i = 0; i < gradient.size(); ++i) {
//     norm += fabs(gradient[i]);
//   }
//   double step = init_step/norm;
//   bool found = false;
//   while (!found && step < param_tol) {
//     step *= 0.5;
//     for (size_t i = 0; i < rates.size() - 1; ++i) {
//       new_rates[i] = rates[i] + gradient[i]*step;
//     }
//     new_rates.back() = new_rates[0]*new_rates[5]*
//       pow(new_rates[3]/ new_rates[1], 2)/new_rates[2];
//     found = true;
//     for (size_t i = 0; i < rates.size(); ++i) {
//       if (new_rates[i] <= 0) found = false;
//     }
//   }
//   abs_diff = step*norm;
// }


// double gradient_ascent(const double param_tol,
//                        const vector<vector<Jump> > &jumps,
//                        const vector<double> &J,
//                        const vector<double> &D,
//                        const vector<double> &rates,
//                        vector<double> &new_rates) {
//   /* compute llk and gradient */
//   vector<double> g;
//   double l = llk(J, D, rates);
//   get_gradient(J, D, rates, g);
//   double step = 0.2;
//   double diff = 1;
//   double new_l = 0;
//   do {
//     step *= 0.5;
//     candidate_rates(param_tol, step, rates, g, new_rates, diff);
//     new_l = llk(J, D, new_rates);
//   } while (new_l <= l && diff > param_tol);
//   return diff;
// }


// double gradient_ascent(const double param_tol,
//                        const vector<Jump> &jumps,
//                        const vector<double> &J,
//                        const vector<double> &D,
//                        const vector<double> &rates,
//                        vector<double> &new_rates) {
//   /* compute llk and gradient */
//   vector<double> g;
//   double l = llk(J, D, rates);
//   get_gradient(J, D, rates, g);
//   double step = 0.2;
//   double diff = 1;
//   double new_l = 0;
//   do {
//     step *= 0.5;
//     candidate_rates(param_tol, step, rates, g, new_rates, diff);
//     new_l = llk(J, D, new_rates);
//   } while (new_l <= l && diff > param_tol);
//   return diff;
// }


// void est_param(const double param_tol,
//                const vector<vector<Jump> > &jumps,
//                const vector<double> &rates,
//                vector<double> &new_rates) {
//   vector<double> J;
//   vector<double> D;
//   get_suff_stat(jumps, J, D);

//   vector<double> current_rates = rates;
//   double current_llk = llk(J, D, current_rates);
//   double diff = gradient_ascent(param_tol, jumps, J, D,
//                                 current_rates, new_rates);
//   double new_llk = llk(J, D, new_rates);

//   while (new_llk - current_llk > param_tol && diff > param_tol) {
//     current_llk = new_llk;
//     current_rates = new_rates;
//     diff = gradient_ascent(param_tol, jumps, J, D,
//                            current_rates, new_rates);
//     new_llk = llk(J, D, new_rates);

//     if (DEBUG) {
//       cerr << "***rates:\t";
//       for (size_t i = 0; i < new_rates.size(); ++i)
//         cerr << new_rates[i] << "\t";
//       cerr << new_llk << endl;
//     }
//   }
// }

// void est_param(const double param_tol,
//                const vector<Jump> &jumps,
//                const vector<double> &rates,
//                vector<double> &new_rates) {
//   vector<double> J;
//   vector<double> D;
//   get_suff_stat(jumps, J, D);

//   vector<double> current_rates = rates;
//   double current_llk = llk(J, D, current_rates);
//   double diff = gradient_ascent(param_tol, jumps, J, D,
//                                 current_rates, new_rates);
//   double new_llk = llk(J, D, new_rates);
//   while (new_llk - current_llk > param_tol && diff > param_tol) {
//     current_llk = new_llk;
//     current_rates = new_rates;
//     diff = gradient_ascent(param_tol, jumps, J, D,
//                            current_rates, new_rates);
//     new_llk = llk(J, D, new_rates);

//     if (DEBUG) {
//       cerr << "***rates:\t";
//       for (size_t i = 0; i < new_rates.size(); ++i)
//         cerr << new_rates[i] << "\t";
//       cerr << new_llk << endl;
//     }
//   }
// }

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
        cerr<< "Branch-" << branch << endl;
        cerr << "Total sites: " << paths[branch].size() << endl;
        cerr << "Total jumps: " << jumps[branch].size() << endl;
      }

      //////////
      if (DEBUG) {
        cerr << "==========" << endl
             << "First 10 jumps:" << endl;
        for (size_t i = 0; i < 10; ++i)
          cerr << "position-" << jumps[branch][i].pos
               << "\ttime-" << jumps[branch][i].t
               << "\tcontext-" << jumps[branch][i].context << endl;
        cerr << "==========" << endl;
      }
      //////////
    }

    if (!param_file.empty()) {
      model_param p;
      read_param(param_file, p);
      vector<double> J;
      vector<double> D;
      get_suff_stat(jumps, J, D);

      //////////
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
      //////////

      vector<double> rates;
      get_rates(p, rates);

      if (VERBOSE) {
        cerr << "starting rates:\t";
        for (size_t i = 0; i < rates.size(); ++i) cerr << rates[i] << "\t";
        cerr << endl;
        cerr << "log-likelihood= " << llk(J, D, rates) << endl;
      }
      if (DEBUG) {
        cerr << "tot_freq:" << endl;
        for (size_t i =0; i < 8; ++i)
          cerr << "[" << i << "]\t" << J[i] << endl;

        cerr << "weights:" << endl;
        for (size_t i =0; i < 8; ++i)
          cerr << "[" << i << "]\t" << D[i] << endl;
      }

      vector<double> gradient;
      get_gradient_logrates(J, D, rates, gradient);

      if (DEBUG) {
        cerr << "gradient :" << endl;
        for (size_t i =0; i < gradient.size(); ++i)
          cerr << "[" << i << "]\t" << gradient[i] << endl;
      }

      double param_tol = 1e-8;
      vector<double> new_rates;
      //est_param(param_tol, jumps, rates, new_rates);
      est_param_logrates(param_tol, jumps, rates, new_rates);

      if (VERBOSE) {
        cerr << "new rates:\t" << endl;
        for (size_t i = 0; i < new_rates.size(); ++i)
          cerr << new_rates[i] << "\t";
        cerr << endl
             <<"new log-likelihood= "
             << llk(J, D, new_rates) << endl;
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
