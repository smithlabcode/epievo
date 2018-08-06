/* Copyright (C) 2018 University of Southern California
 *                    Jianghan Qu, Andrew D Smith and Xiaojing Ji
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
#ifndef END_COND_SAMPLING_HPP
#define END_COND_SAMPLING_HPP

#include "ContinuousTimeMarkovModel.hpp"
#include "Path.hpp"
#include <vector>

void
end_cond_sample_direct(const CTMarkovModel &the_model,
                       const size_t start_state, const size_t end_state,
                       const double time_interval,
                       std::mt19937 &gen, std::vector<double> &jump_times,
                       const double start_time = 0.0);

bool
end_cond_sample_forward_rejection(const size_t max_sample_count,
                                  const CTMarkovModel &the_model,
                                  const size_t start_state, const size_t end_state,
                                  const double time_interval,
                                  std::mt19937 &gen, std::vector<double> &jump_times,
                                  const double start_time = 0.0);

void
end_cond_sample_unif(const CTMarkovModel &the_model,
                     const size_t start_state, const size_t end_state,
                     const double time_interval,
                     std::mt19937 &gen, std::vector<double> &jump_times,
                     std::vector<mixJump> &mjumps,
                     const double start_time = 0.0);


bool
end_cond_sampling_Nielsen(const size_t max_sample_count,
                          const CTMarkovModel &the_model,
                          const size_t start_state, const size_t end_state,
                          const double time_interval,
                          std::mt19937 &gen, std::vector<double> &jump_times,
                          const double start_time = 0.0);

double
end_cond_sample_prob(const CTMarkovModel &the_model,
                     const std::vector<double> &jump_times,
                     const size_t start_state, const size_t end_state,
                     const double start_time, const double end_time,
                     size_t start_jump, const size_t end_jump);

double
end_cond_sample_unif_prob(const CTMarkovModel &the_model,
                          const std::vector<mixJump> &mjumps,
                          const size_t start_state, const size_t end_state,
                          const double start_time, const double end_time,
                          size_t start_jump, const size_t end_jump);

double
forward_sample_prob(const CTMarkovModel &the_model,
                     const std::vector<double> &jump_times,
                     const size_t start_state, const size_t end_state,
                     const double start_time, const double end_time,
                     size_t start_jump, const size_t end_jump);

#endif
