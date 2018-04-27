#include "TripletPattern.hpp"

#include <vector>
#include <random>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cassert>
#include <numeric>  /* std::inner_product */


using std::vector;
using std::cerr;
using std::endl;

/* PatSeq */
PatSeq::PatSeq(const vector<bool> &seq) {
  const size_t n = seq.size();
  vector<vector<size_t> > pos_grp(8, vector<size_t>() );
  for (size_t i = 1; i < n-1; ++i) {
    size_t pat = (size_t)(seq[i-1])*4 + (size_t)(seq[i])*2 + (size_t)(seq[i+1]);
    pos_grp[pat].push_back(i);
  }

  /* set pat_freq and pos_by_pat */
  pos_by_pat.resize(2, 0);
  cum_pat_freq.resize(8, 0);
  for (size_t i = 0; i < 8; ++i) {
    pos_by_pat.insert(pos_by_pat.end()-1, pos_grp[i].begin(), pos_grp[i].end());
    if (i == 0) {
      cum_pat_freq[i] = pos_grp[i].size();
    } else {
      cum_pat_freq[i] = cum_pat_freq[i-1] + pos_grp[i].size();
    }
  }

  /* set idx_in_pat */
  idx_in_pat.resize(n, 0);
  for (size_t i = 0; i < n; ++i) {
    if (i==0 || i==n-1) {
      idx_in_pat[i] = i;
    } else {
      size_t pos = pos_by_pat[i];
      idx_in_pat[pos] = i;
    }
  }

  start_state = seq[0];
  end_state = seq.back();
}

void PatSeq::get_all_context_freq(vector<size_t> &freq) const {
  freq.resize(8, 0);
  freq[0] = cum_pat_freq[0];
  for (size_t i = 1; i < 8; ++i)
    freq[i] = cum_pat_freq[i] - cum_pat_freq[i-1];
}


size_t PatSeq::get_context_freq(const size_t context) const {
  assert(context < 8);
  if (context == 0)
    return cum_pat_freq[0];
  else
    return cum_pat_freq[context] - cum_pat_freq[context-1];
}


size_t PatSeq::get_context(const size_t pos) const {
  if (pos == 0 || pos == pos_by_pat.size() - 1) {
    return 8;
  } else {
    size_t idx = idx_in_pat[pos];
    size_t pat = 0;
    while (cum_pat_freq[pat] < idx) ++pat;
    return pat;
  }
}

void PatSeq::single_update(const size_t pos, const size_t context,
                           const size_t to_context) {
  size_t loc = idx_in_pat[pos];
  if (to_context < context) {
    const size_t block_start = cum_pat_freq[context -1] + 1;
    // swap to beginning of chunck
    iter_swap(pos_by_pat.begin() + block_start, pos_by_pat.begin() + loc);
    iter_swap(idx_in_pat.begin() + pos_by_pat[block_start],
              idx_in_pat.begin() + pos_by_pat[loc]);
    size_t cur_context = context;
    while (cur_context > to_context) {
      size_t prev_context = cur_context - 1;
      size_t cur_block_start = cum_pat_freq[prev_context] + 1;
      size_t prev_block_start = (prev_context > 0)? (cum_pat_freq[prev_context-1] + 1) : 1;
      iter_swap(pos_by_pat.begin() + cur_block_start,
                pos_by_pat.begin() + prev_block_start);
      iter_swap(idx_in_pat.begin() + pos_by_pat[cur_block_start],
                idx_in_pat.begin() + pos_by_pat[prev_block_start]);

      /* update cum_pat_freq */
      ++cum_pat_freq[prev_context];
      --cur_context;
    }
  } else {
    /* swap to end of chunck */
    const size_t block_end = cum_pat_freq[context];
    iter_swap(pos_by_pat.begin() + block_end, pos_by_pat.begin() + loc);
    iter_swap(idx_in_pat.begin() + pos_by_pat[block_end],
              idx_in_pat.begin() + pos_by_pat[loc]);
    size_t cur_context = context;
    while (cur_context < to_context) {
      size_t next_context = cur_context + 1;
      size_t cur_block_end = cum_pat_freq[cur_context];
      size_t next_block_end = cum_pat_freq[next_context];
      iter_swap(pos_by_pat.begin() + cur_block_end,
                pos_by_pat.begin() + next_block_end);
      iter_swap(idx_in_pat.begin() + pos_by_pat[cur_block_end],
                idx_in_pat.begin() + pos_by_pat[next_block_end]);
      /* update cum_pat_freq */
      --cum_pat_freq[cur_context];
      ++cur_context;
    }
  }
}

size_t PatSeq::random_mutate(const size_t context, std::mt19937 &gen) {
  assert(context < 8);
  size_t loc = (context ==0) ? 0 : cum_pat_freq[context - 1];
  size_t pat_freq =
    (context == 0)? cum_pat_freq[context] : cum_pat_freq[context] - cum_pat_freq[context - 1];

  std::uniform_int_distribution<int> dist(1, pat_freq);
  loc += dist(gen);

  const size_t pos = pos_by_pat[loc];
  const size_t pos_l = pos - 1;
  const size_t context_l = get_context(pos_l);
  const size_t pos_r = pos + 1;
  const size_t context_r = get_context(pos_r);

  single_update(pos, context, context^2);  // flip the middle bit
  if (pos_l > 0)
    single_update(pos_l, context_l, context_l^1);  // flip the right bit
  if (pos_r < pos_by_pat.size() - 1)
    single_update(pos_r, context_r, context_r^4);  // flip the left bit

  return pos;
}


void PatSeq::mutate(const size_t pos) {
  const size_t context = get_context(pos);
  const size_t pos_l = pos - 1;
  const size_t context_l = get_context(pos_l);
  const size_t pos_r = pos + 1;
  const size_t context_r = get_context(pos_r);

  single_update(pos, context, context^2);  // flip the middle bit
  if (pos_l > 0)
    single_update(pos_l, context_l, context_l^1);  // flip the right bit
  if (pos_r < pos_by_pat.size() - 1)
    single_update(pos_r, context_r, context_r^4);  // flip the left bit
}


void PatSeq::to_seq(std::vector<bool> &seq) const {
  size_t n_site = idx_in_pat.size();
  seq.resize(n_site, true);
  size_t idx = 1;
  for (size_t pat = 0; pat < 8; ++pat) {
    bool state = (pat >>1) % 2;
    while (idx <= cum_pat_freq[pat]) {
      seq[pos_by_pat[idx]] = state;
      ++idx;
    }
  }
  seq[0] = start_state;
  seq.back() = end_state;
}


void PatSeq::to_domain_sizes(std::vector<size_t> &domain_sizes) const {
  domain_sizes.clear();
  std::vector<bool> seq;
  to_seq(seq);
  size_t s = 0;
  bool in_domain = false;
  for (size_t i = 0; i < seq.size(); ++i) {
    if (seq[i] && !in_domain) {
      s = 1;
      in_domain = true;
    } else if (seq[i] && in_domain) {
      ++s;
    } else if (!seq[i] && in_domain) {
      in_domain = false;
      domain_sizes.push_back(s);
      s = 0;
    }
  }
}

