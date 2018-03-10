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

#include "TripletSampler.hpp"

#include "StateSeq.hpp"

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>
#include <cassert>

using std::vector;

TripletSampler::TripletSampler(const vector<char> &seq) {
  assert(seq.size() >= 3); // must have at least one triplet

  const size_t seq_len = seq.size();

  // get the cumulative counts of each binary triplet
  cum_pat_count.resize(8, 0ul);
  for (size_t i = 1; i < seq_len - 1; ++i)
    ++cum_pat_count[triple2idx(seq[i-1], seq[i], seq[i+1])];
  std::partial_sum(cum_pat_count.begin(), cum_pat_count.end(),
                   cum_pat_count.begin());

  idx_in_pat = vector<size_t>(seq_len);  // valid for all sequence positions
  pos_by_pat = vector<size_t>(seq_len - 2); // total number of triples

  // place each sequence position into the appropriate part of the
  // vector; this is just one pass of counting-sort
  for (size_t i = 1; i < seq_len - 1; ++i) {
    const size_t idx = triple2idx(seq[i-1], seq[i], seq[i+1]);
    --cum_pat_count[idx];
    const size_t x = cum_pat_count[idx];
    pos_by_pat[x] = i;
    idx_in_pat[i] = x;
  }

  /* notice that cum_pat_count has 9 elements */
  cum_pat_count.push_back(seq_len - 2); // total number of triples

  // remember the start and end state; these will never change
  start_state = seq[0];
  end_state = seq.back();
}


TripletSampler::TripletSampler(const StateSeq &s) {
  assert(s.seq.size() >= 3); // must have at least one triplet

  const size_t seq_len = s.seq.size();

  // get the cumulative counts of each binary triplet
  cum_pat_count.resize(8, 0ul);
  for (size_t i = 1; i < seq_len - 1; ++i)
    ++cum_pat_count[triple2idx(s.seq[i-1], s.seq[i], s.seq[i+1])];
  std::partial_sum(cum_pat_count.begin(), cum_pat_count.end(),
                   cum_pat_count.begin());

  idx_in_pat = vector<size_t>(seq_len);  // valid for all sequence positions
  pos_by_pat = vector<size_t>(seq_len - 2); // total number of triples

  // place each sequence position into the appropriate part of the
  // vector; this is just one pass of counting-sort
  for (size_t i = 1; i < seq_len - 1; ++i) {
    const size_t idx = triple2idx(s.seq[i-1], s.seq[i], s.seq[i+1]);
    --cum_pat_count[idx];
    const size_t x = cum_pat_count[idx];
    pos_by_pat[x] = i;
    idx_in_pat[i] = x;
  }

  /* notice that cum_pat_count has 9 elements */
  cum_pat_count.push_back(seq_len - 2); // total number of triples

  // remember the start and end state; these will never change
  start_state = s.seq[0];
  end_state = s.seq.back();
}


void
TripletSampler::get_triplet_counts(vector<size_t> &triplet_counts) const {
  triplet_counts.clear();
  adjacent_difference(cum_pat_count.begin() + 1, cum_pat_count.end(),
                      std::back_inserter(triplet_counts));
}

size_t
TripletSampler::get_triplet_count(const size_t context) const {
  assert(context < 8);
  return cum_pat_count[context + 1] - cum_pat_count[context];
}


void
TripletSampler::single_update(const size_t pos, size_t context,
                              const size_t to_context) {

  size_t loc = idx_in_pat[pos];
  if (context > to_context) {

    size_t block_start = cum_pat_count[context];

    // swap current location to beginning of block
    iter_swap(pos_by_pat.begin() + block_start, pos_by_pat.begin() + loc);
    iter_swap(idx_in_pat.begin() + pos_by_pat[block_start],
              idx_in_pat.begin() + pos_by_pat[loc]);

    // swap beginning elements of blocks
    while (context > to_context + 1) {
      ++cum_pat_count[context]; // expand block for preceding context
      --context; // move to preceding context

      /* swap last to first within block for preceding context */
      const size_t preceding_block_start = cum_pat_count[context];
      iter_swap(pos_by_pat.begin() + block_start,
                pos_by_pat.begin() + preceding_block_start);
      iter_swap(idx_in_pat.begin() + pos_by_pat[block_start],
                idx_in_pat.begin() + pos_by_pat[preceding_block_start]);
      block_start = preceding_block_start;
    }
    ++cum_pat_count[context]; // no swap for last context (block boundary moves)
  }
  else {

    /* (ADS) for this case, I originally coded it using a "context +
       1" in most places, but then realized it would be the same if I
       just incremented the context at the start: */
    ++context;

    size_t block_end = cum_pat_count[context] - 1; // last position in
                                                   // current block

    // swap current location to end of block
    iter_swap(pos_by_pat.begin() + block_end, pos_by_pat.begin() + loc);
    iter_swap(idx_in_pat.begin() + pos_by_pat[block_end],
              idx_in_pat.begin() + pos_by_pat[loc]);

    // swap ending elements of blocks
    while (context < to_context) {
      --cum_pat_count[context]; // grow new context by shifting boundary
      ++context; // now move to next block

      /* swap first-to-last within block for next context */
      const size_t next_block_end = cum_pat_count[context] - 1;
      iter_swap(pos_by_pat.begin() + block_end,
                pos_by_pat.begin() + next_block_end);
      iter_swap(idx_in_pat.begin() + pos_by_pat[block_end],
                idx_in_pat.begin() + pos_by_pat[next_block_end]);
      block_end = next_block_end;
    }
    --cum_pat_count[context];
  }
}


/* get_context is used to lookup the context of a particular triplet
   using its position in the pos_by_pat. The vector pos_by_pat is
   organized such that all positions having a particular pattern
   appear consecutively, but are otherwise unordered.
*/
size_t
TripletSampler::get_context(const size_t pos) const {
  assert(pos > 0 && pos < idx_in_pat.size() - 1);
  const size_t idx = idx_in_pat[pos];
  size_t pat = 0;
  // below: iterates at most 8 times; using "pat + 1" because the
  // cum_pat_count is shifted. Also, the "<=" is needed in case some
  // set of triplets immediately preceding the triplet at pos appear 0
  // times. For example, if there are no 000 or 001 patterns, and pos
  // appears first among the 010 patterns, then cum_pat_count[pat + 1]
  // would be equal to idx immediately, and we would return
  // pat=0. This can only happen in degenerate cases.
  // JQU: cum_pat_count[pat + 1] - 1 is the index of the last element
  // in the block of pat in pos_by_pat, so advance pat by 1 as long
  // as cum_pat_count[pat + 1] - 1 < idx
  while (cum_pat_count[pat + 1] <= idx) ++pat;
  return pat;
}


void
TripletSampler::mutate(const size_t pos, const size_t context) {

  // flip the middle bit of the current triplet pattern
  single_update(pos, context, flip_mid_bit(context));

  // flip the right bit of the left neighbor's triplet pattern (unless
  // the left neighbor is not part of a full triplet due to being at
  // the start of the sequence)
  const size_t pos_l = pos - 1;
  if (pos_l > 0) {
    const size_t context_l = get_context(pos_l);
    single_update(pos_l, context_l, flip_right_bit(context_l));
  }

  // flip the left bit of the right neighbor's triplet pattern (unless
  // the right neighbor would be past the end of the sequence).
  const size_t pos_r = pos + 1;
  if (pos_r < pos_by_pat.size() - 1) {
    const size_t context_r = get_context(pos_r);
    single_update(pos_r, context_r, flip_left_bit(context_r));
  }
}


size_t
TripletSampler::random_mutate(const size_t context, std::mt19937 &gen) {
  assert(context < 8);

  // below subtracting 1 because sampling is closed interval
  std::uniform_int_distribution<size_t> dist(cum_pat_count[context],
                                             cum_pat_count[context + 1] - 1);
  const size_t pos = pos_by_pat[dist(gen)];
  mutate(pos, context);

  return pos;
}


// extract the sequence by "unpermuting" the positions
void
TripletSampler::get_sequence(std::vector<char> &seq) const {
  seq.resize(idx_in_pat.size(), true);

  size_t idx = 1;
  for (size_t pat = 0; pat < 8; ++pat) {
    const bool state = get_mid_bit(pat);
    while (idx <= cum_pat_count[pat]) {
      seq[pos_by_pat[idx]] = state;
      ++idx;
    }
  }

  seq[0] = start_state;
  seq.back() = end_state;
}

// same as above, but for the StateSeq object
void
TripletSampler::get_sequence(StateSeq &s) const {
  get_sequence(s.seq);
}
