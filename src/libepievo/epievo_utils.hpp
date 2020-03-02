/* Copyright (C) 2020 University of Southern California
 *                    Xiaojing Ji and Andrew D Smith
 *
 * Author: Andrew D. Smith and Xiaojing Ji
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

#ifndef EPIEVO_UTILS_HPP
#define EPIEVO_UTILS_HPP

#include <vector>
#include <string>
#include <cassert>

struct two_by_two {

  // constructor
  two_by_two() : v00(0.0), v01(0.0), v10(0.0), v11(0.0) {};
  two_by_two(const double w00, const double w01,
             const double w10, const double w11) :
    v00(w00), v01(w01), v10(w10), v11(w11) {};

  // accessing an element
  constexpr double operator() (const size_t r, const size_t c) const {
    return (r == 0 ? (c == 0 ? v00 : v01) : (c == 0 ? v10 : v11));
  }
  // changing an element
  double &operator() (const size_t r, const size_t c) {
    return (r == 0 ? (c == 0 ? v00 : v01) : (c == 0 ? v10 : v11));
  }
  // get a row
  std::vector<double> operator[](const size_t r) const {
    return r == 0 ? std::vector<double>{v00, v01} :
    std::vector<double>{v10, v11};
  }
  two_by_two &operator+=(const two_by_two &rhs) {
    v00 += rhs.v00; v01 += rhs.v01;
    v10 += rhs.v10; v11 += rhs.v11;
    return *this
  }
  void divide_all(const double &denom) {
    v00 /= denom; v01 /= denom;
    v10 /= denom; v11 /= denom;
  }
  // reset all values to zero
  void reset() {v00 = 0.0; v01 = 0.0; v10 = 0.0; v11 = 0.0;}
  double v00, v01, v10, v11;
};

typedef std::vector<bool> state_seq;

void
get_triplet_counts(const state_seq &s,
                   std::vector<size_t> &triplet_counts);
void
get_triplet_proportions(const state_seq &s,
                        std::vector<double> &triplet_proportions);
void
get_pair_counts(const state_seq &s, std::vector<size_t> &pair_counts);
void
get_pair_proportions(const state_seq &s, std::vector<double> &pair_proportions);

std::string
summary_string(const state_seq &s);

void
read_states_file(const std::string &statesfile, std::vector<std::string> &names,
                 std::vector<state_seq> &the_states);

// utilities for working with binary states, especially triplets
constexpr size_t
triple2idx(const bool i, const bool j, const bool k) {
  return i*4 + j*2 + k;
}

constexpr size_t
pair2idx(const bool i, const bool j) {
  return i*2 + j;
}

// complementing a state
constexpr size_t
complement_state(const size_t x) {
  // assert(x == 0ul || x == 1ul);
  return 1ul - x;
}

constexpr bool
complement_state(const bool x) {
  return !x;
}

// flipping bits
constexpr size_t
flip_left_bit(const size_t x)  {return x ^ 4ul;}
constexpr size_t
flip_mid_bit(const size_t x)   {return x ^ 2ul;}
constexpr size_t
flip_right_bit(const size_t x) {return x ^ 1ul;}

// accessing bits
constexpr bool get_left_bit(const size_t x)  {return x & 4ul;}
constexpr bool get_mid_bit(const size_t x)   {return x & 2ul;}
constexpr bool get_right_bit(const size_t x) {return x & 1ul;}

constexpr bool get_left_bit_from_pair(const size_t x)  {return x & 2ul;}
constexpr bool get_right_bit_from_pair(const size_t x) {return x & 1ul;}

inline void
get_bits_from_triple(const size_t x, bool &l, bool &m, bool &r) {
  l = get_left_bit(x);
  m = get_mid_bit(x);
  r = get_right_bit(x);
}

inline void
get_bits_from_pair(const size_t x, bool &l, bool &r) {
  l = get_left_bit_from_pair(x);
  r = get_right_bit_from_pair(x);
}

#include <bitset>
#include <sstream>

template <class T>
std::string
triplet_info_to_string(const std::vector<T> &v) {
  static const size_t n_triplets = 8;
  assert(v.size() >= n_triplets);
  std::ostringstream oss;
  oss << std::bitset<3>(0) << '\t' << v.front();
  for (size_t i = 1; i < n_triplets; ++i)
    oss << '\n' << std::bitset<3>(i) << '\t' << v[i];
  return oss.str();
}

template <class T>
std::string
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
