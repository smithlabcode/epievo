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

#include "StateSeq.hpp"

#include <vector>
#include <cassert>
#include <functional>
#include <algorithm>

using std::vector;

void
StateSeq::get_domain_sizes(vector<size_t> &domain_sizes) const {
  domain_sizes.clear();
  size_t s = 0;
  bool in_domain = false;
  for (size_t i = 0; i < seq.size(); ++i) {
    if (seq[i] && !in_domain) {
      s = 1;
      in_domain = true;
    }
    else if (seq[i] && in_domain) {
      ++s;
    }
    else if (!seq[i] && in_domain) {
      in_domain = false;
      domain_sizes.push_back(s);
      s = 0;
    }
  }
}

void
StateSeq::get_triplet_counts(std::vector<size_t> &triplet_counts) const {
  static const size_t n_triplets = 8;

  triplet_counts.resize(n_triplets, 0);
  for (size_t i = 2; i < seq.size(); ++i)
    triplet_counts[triple2idx(seq[i-2], seq[i-1], seq[i])]++;
}

void
StateSeq::get_triplet_proportions(std::vector<double> &triplet_props) const {

  vector<size_t> triplet_counts;
  get_triplet_counts(triplet_counts);

  triplet_props.resize(triplet_counts.size());
  std::transform(triplet_counts.begin(), triplet_counts.end(),
                 triplet_props.begin(), std::bind2nd(std::divides<double>(),
                                                     seq.size() - 2));
}
