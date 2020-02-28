/* Copyright (C) 2019 University of Southern California
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

#ifndef TRIPLET_SAMPLER_HPP
#define TRIPLET_SAMPLER_HPP

#include <vector>
#include <random>
#include <algorithm>

#include "epievo_utils.hpp"

class TripletSampler {
public:
  TripletSampler(const state_seq &seq);
  void get_triplet_counts(std::vector<size_t> &counts) const;
  size_t get_triplet_count(const size_t context) const;
  size_t random_mutate(const size_t context, std::mt19937 &gen);

  // mutate: the second argument can be computed using the first, but
  // passing it as parameter skips having to compute it again
  void mutate(const size_t pos, const size_t context);
  void mutate(const size_t pos) {mutate(pos, get_context(pos));}

  // get_sequence: extracts the sequences corresponding to the current
  // state of the calling TripletSampler object, which is obtained by
  // "unpermuting" the idx_in_pat; essentially placing the correct bit
  // into each position of the sequence.
  void get_sequence(state_seq &s) const;

private:
  /* positions in state sequence organized by triplet pattern */
  std::vector<size_t> pos_by_pat;

  /* indices in pos_by_pat for positions in binary-state sequence*/
  std::vector<size_t> idx_in_pat;

  /* number of occurrences of triplets preceding current one; has 9
     elements and the last one is the total number of all triplets */
  std::vector<size_t> cum_pat_count;

  bool start_state;
  bool end_state;

  size_t get_context(const size_t pos) const;
  void single_update(const size_t pos, const size_t context,
                     const size_t to_context);
};

#endif
