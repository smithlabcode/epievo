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

struct StateSeq {
  StateSeq(const std::vector<char> &s) : seq(s) {}
  void to_domain_sizes(std::vector<size_t> &domain_sizes) const;
  std::vector<char> seq;
};

template <class T> size_t
triple2idx(const T i, const T j, const T k) {
  return i*4 + j*2 + k;
}

inline size_t
flip_right_bit(const size_t x) {return x ^ 1ul;}

inline size_t
flip_mid_bit(const size_t x) {return x ^ 2ul;}

inline size_t
flip_left_bit(const size_t x) {return x ^ 4ul;}

inline bool
get_right_bit(const size_t x) {return x & 1ul;}

inline bool
get_mid_bit(const size_t x) {return x & 2ul;}

inline bool
get_left_bit(const size_t x) {return x & 4ul;}

#endif
