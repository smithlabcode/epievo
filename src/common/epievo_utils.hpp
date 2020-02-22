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
  // reset all values to zero
  void reset() {v00 = 0.0; v01 = 0.0; v10 = 0.0; v11 = 0.0;}
  double v00, v01, v10, v11;
};

#endif
