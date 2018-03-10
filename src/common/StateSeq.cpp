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
#include <string>
#include <bitset>
#include <sstream>

using std::vector;
using std::string;

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

void
StateSeq::get_pair_counts(std::vector<size_t> &pair_counts) const {
  static const size_t n_pairs = 4;

  pair_counts.resize(n_pairs, 0);
  for (size_t i = 1; i < seq.size(); ++i)
    pair_counts[pair2idx(seq[i-1], seq[i])]++;
}

void
StateSeq::get_pair_proportions(std::vector<double> &pair_props) const {

  vector<size_t> pair_counts;
  get_pair_counts(pair_counts);

  pair_props.resize(pair_counts.size());
  std::transform(pair_counts.begin(), pair_counts.end(),
                 pair_props.begin(), std::bind2nd(std::divides<double>(),
                                                  seq.size() - 1));
}


template <class T> static string
triplet_info_to_string(const std::vector<T> &v) {
  static const size_t n_triplets = 8;
  assert(v.size() >= n_triplets);
  std::ostringstream oss;
  oss << std::bitset<3>(0) << '\t' << v.front();
  for (size_t i = 1; i < n_triplets; ++i)
    oss << '\n' << std::bitset<3>(i) << '\t' << v[i];
  return oss.str();
}


template <class T> static string
pair_info_to_string(const std::vector<T> &v) {
  static const size_t n_pairs = 4;

  assert(v.size() >= n_pairs);
  std::ostringstream oss;
  oss << std::bitset<2>(0) << '\t' << v.front();
  for (size_t i = 1; i < n_pairs; ++i)
    oss << '\n' << std::bitset<2>(i) << '\t' << v[i];
  return oss.str();
}

string
StateSeq::summary_string() const {

  vector<double> triplet_props;
  get_triplet_proportions(triplet_props);
  vector<double> pair_props;
  get_pair_proportions(pair_props);

  std::ostringstream oss;
  oss << "triplet proportions:\n"
      << triplet_info_to_string(triplet_props) << '\n'
      << "pair proportions:\n"
      << pair_info_to_string(pair_props);
  return oss.str();
}
