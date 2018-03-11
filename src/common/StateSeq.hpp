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

#ifndef STATESEQ_HPP
#define STATESEQ_HPP

#include <vector>
#include <cstdlib>
#include <string>

struct StateSeq {

  StateSeq() {}
  StateSeq(const std::vector<char> &s) : seq(s) {}
  void get_domain_sizes(std::vector<size_t> &domain_sizes) const;
  std::vector<char> seq;

  void get_triplet_counts(std::vector<size_t> &triplet_counts) const;
  void get_triplet_proportions(std::vector<double> &triplet_proportions) const;

  void get_pair_counts(std::vector<size_t> &pair_counts) const;
  void get_pair_proportions(std::vector<double> &pair_proportions) const;

  std::string summary_string() const;
};

inline size_t
triple2idx(const bool i, const bool j, const bool k) {
  return i*4 + j*2 + k;
}

inline size_t
pair2idx(const bool i, const bool j) {
  return i*2 + j;
}

// flipping bits
inline size_t flip_left_bit(const size_t x)  {return x ^ 4ul;}
inline size_t flip_mid_bit(const size_t x)   {return x ^ 2ul;}
inline size_t flip_right_bit(const size_t x) {return x ^ 1ul;}

// accessing bits
inline bool get_left_bit(const size_t x)  {return x & 4ul;}
inline bool get_mid_bit(const size_t x)   {return x & 2ul;}
inline bool get_right_bit(const size_t x) {return x & 1ul;}

inline void
get_bits_from_triple(const size_t x, bool &l, bool &m, bool &r) {
  l = get_left_bit(x);
  m = get_mid_bit(x);
  r = get_right_bit(x);
}

#include <bitset>
#include <sstream>
#include <cassert>

template <class T> std::string
triplet_info_to_string(const std::vector<T> &v) {
  static const size_t n_triplets = 8;
  assert(v.size() >= n_triplets);
  std::ostringstream oss;
  oss << std::bitset<3>(0) << '\t' << v.front();
  for (size_t i = 1; i < n_triplets; ++i)
    oss << '\n' << std::bitset<3>(i) << '\t' << v[i];
  return oss.str();
}


template <class T> std::string
pair_info_to_string(const std::vector<T> &v) {
  static const size_t n_pairs = 4;

  assert(v.size() >= n_pairs);
  std::ostringstream oss;
  oss << std::bitset<2>(0) << '\t' << v.front();
  for (size_t i = 1; i < n_pairs; ++i)
    oss << '\n' << std::bitset<2>(i) << '\t' << v[i];
  return oss.str();
}

#endif
