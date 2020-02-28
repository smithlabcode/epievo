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
 */

#ifndef EMIT_DISTRO_HPP
#define EMIT_DISTRO_HPP

#include <vector>

struct Bernoulli {
  Bernoulli() : p(0.5) {}
  Bernoulli(const double p_) : p(p_) {}
  double operator()(const bool val) const {
    return val ? p : 1.0 - p;
  }
  void fit(const std::vector<bool> &vals) {
    p = std::accumulate(std::begin(vals), std::end(vals), 0.0)/vals.size();
  }
  double p;
};

#endif
