/* Copyright (C) 2020 University of Southern California
 *                    Xiaojing Ji, Jianghan Qu and Andrew D Smith
 *
 * Author: Andrew D. Smith, Jianghan Qu and Xiaojing Ji
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

#include <string>
#include <vector>
#include <cassert>

using std::vector;

template <typename T = double>
struct two_by_two {
  
  // constructor
  two_by_two() : v00(0.0), v01(0.0), v10(0.0), v11(0.0) {};
  two_by_two(const T w00, const T w01, const T w10, const T w11) :
    v00(w00), v01(w01), v10(w10), v11(w11) {};

  // accessing an element
  constexpr T operator() (const size_t r, const size_t c) const {
    return (r == 0 ? (c == 0 ? v00 : v01) : (c == 0 ? v10 : v11));
  };
  
  // changing an element
  T &operator() (const size_t r, const size_t c) {
    return (r == 0 ? (c == 0 ? v00 : v01) : (c == 0 ? v10 : v11));
  };
  
  // get a row
  vector<T> operator[] (const size_t r) const {
    return (r == 0? vector<double> {v00, v01} :
                    vector<double> {v10, v11});
  }
  
  // reset all values to zero
  void reset() {
    v00 = 0.0;
    v01 = 0.0;
    v10=0.0;
    v11=0.0;
  }

  T v00;
  T v01;
  T v10;
  T v11;
};

#endif
