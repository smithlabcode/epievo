
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <math.h>   /* exp, sqrt, pow */
#include <numeric>  /* std::inner_product */
#include <algorithm>    // std::sort

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "PhyloTreePreorder.hpp"
#include "path.hpp"  /* related to path and vector of paths*/
#include "TwoStateSeq.hpp"

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
  Jump(double time, size_t position): t(time), pos(position) {};
};

bool jump_compare(Jump lhs, Jump rhs) { return lhs.t < rhs.t;}


// full paths -> initial seq, jumps -> sorted jumps -> jump context

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
  TwoStateSeq tss(seq);
  vector<vector<size_t> > freq_before_jump;

  for (size_t i =0; i < jumps.size(); ++i) {
    vector<size_t> freq;
    tss.count_s3(freq);
    freq_before_jump.push_back(freq);
    tss.mutate(jumps[i].pos, jumps[i].context); //get mutation context
  }
}

void est_two_state_markov_trans_prob(const vector<bool> &seq,
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
    bool VERBOSE = false;


    OptionParser opt_parse(strip_path(argv[0]), "estimate parameter from paths",
                           "<path-file>");
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

    if (DEBUG) {
    cerr << paths.size() << endl
         << paths[0].size() << "\t positions" << endl
         << paths[0][1].jumps.size() << "\t jumps" << endl;
    }

    const size_t n_sites = paths[0].size();
    vector<bool> seq;
    get_inital_seq(paths[0], seq);

    vector<vector<double> > trans_prob;
    est_two_state_markov_trans_prob(seq, trans_prob);

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

    vector<Jump> jumps;
    get_jumps(paths[0], jumps);
    if (DEBUG) {
      cerr << "[DEBUG]===" << endl
           << "Total jumps: " << jumps.size() << endl
           << "First 10 jumps:" << endl;
      for (size_t i = 0; i < 10; ++i)
        cerr << "position-" << jumps[i].pos << "\ttime-" << jumps[i].t
             << "\tcontext-" << jumps[i].context << endl;
      cerr << "==========" << endl;
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
