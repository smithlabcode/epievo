#include "jump.hpp"

#include <algorithm>    // std::sort
#include <iostream>

using std::vector;

bool
jump_compare(Jump lhs, Jump rhs) { return lhs.t < rhs.t;}

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

void get_scaled_jumps(const vector<Path> &paths, vector<Jump> &jumps) {
  /*only set jumping times and positions */
  for (size_t i = 0; i < paths.size(); ++i) {
    const double tot_time = paths[i].tot_time;
    for (size_t j = 0; j < paths[i].jumps.size(); ++j) {
      jumps.push_back(Jump(paths[i].jumps[j]/tot_time, i));
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


////////////////////////////////////////////////////////////////////////////////
// Estimate rates given branch lenghts and jumping times
////////////////////////////////////////////////////////////////////////////////

void get_suff_stat(const vector<vector<Jump> > &jumps,
                   vector<double> &J,
                   vector<double> &D) {
  /* J_{ijk} = Total number of jumps in context ijk */
  /* D_{ijk} = Sum over all jumps:
     (holding time before such a jump)*(frequency of jump context) */
  J.resize(8, 0);
  D.resize(8, 0);
  for (size_t b = 0; b < jumps.size(); ++b) {
    double prev_jump_time = 0;
    const size_t n_jumps = jumps[b].size();
    for (size_t i = 0; i < n_jumps; ++i) { /*over all jumps*/
      const size_t context = jumps[b][i].context;
      const vector<size_t> freq = jumps[b][i].freq;
      const double jump_time = jumps[b][i].t;
      const double holding_time = jump_time - prev_jump_time;
      J[context] += 1.0;
      for (size_t ct = 0; ct < 8; ++ct) { /*over all contexts*/
        D[ct] += freq[ct]*holding_time;
      }
      prev_jump_time = jump_time;
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
// Estimate rates and branch length together
////////////////////////////////////////////////////////////////////////////////

void get_suff_stat_by_branch(const vector<vector<Jump> > &jumps,
                             vector<double> &J, vector<vector<double> > &D) {
  /* */
  /* J_{ijk} = Total number of jumps in context ijk*/
  /* D_{ijk}{b} = Sum over all jumps on branch b:
     (holding time before such a jump)*(frequency of jump context) */
  J.resize(8, 0);
  D.resize(jumps.size(), vector<double>(8, 0));
  for (size_t b = 0; b < jumps.size(); ++b) {
    double prev_jump_time = 0;
    const size_t n_jumps = jumps[b].size();
    for (size_t i = 0; i < n_jumps; ++i) { /*over all jumps*/
      const size_t context = jumps[b][i].context;
      const vector<size_t> freq = jumps[b][i].freq;
      const double jump_time = jumps[b][i].t;
      const double holding_time = jump_time - prev_jump_time;
      J[context] += 1.0;
      for (size_t ct = 0; ct < 8; ++ct) { /*over all contexts*/
        D[b][ct] += freq[ct]*holding_time;
      }
      prev_jump_time = jump_time;
    }
  }
}
