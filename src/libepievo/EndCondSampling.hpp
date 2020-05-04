/* Copyright (C) 2019 University of Southern California
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

#ifndef END_COND_SAMPLING_HPP
#define END_COND_SAMPLING_HPP

#include "ContinuousTimeMarkovModel.hpp"
#include "Path.hpp"
#include <vector>

/* "direct" */
void
end_cond_sample_direct(const TwoStateCTMarkovModel &the_model,
                       const size_t start_state, const size_t end_state,
                       const double time_interval,
                       std::mt19937 &gen, std::vector<double> &jump_times,
                       const double start_time = 0.0);

/* "forward" */
bool
end_cond_sample_forward_rejection(std::exponential_distribution<double> &exp1,
                                  std::exponential_distribution<double> &exp2,
                                  const size_t start_state, const size_t end_state,
                                  const double time_interval,
                                  std::mt19937 &gen, std::vector<double> &jump_times,
                                  const double start_time = 0.0);

bool
end_cond_sample_forward_rejection(const TwoStateCTMarkovModel &the_model,
                                  const size_t start_state, const size_t end_state,
                                  const double time_interval,
                                  std::mt19937 &gen, std::vector<double> &jump_times,
                                  const double start_time = 0.0);

bool
end_cond_sample_forward_rejection(const TwoStateCTMarkovModel &the_model,
                                  const size_t start_state, const size_t end_state,
                                  const double time_interval,
                                  std::mt19937 &gen,
                                  std::uniform_real_distribution<double> &unif,
                                  std::vector<double> &jump_times,
                                  const double start_time = 0.0);

/* "nielsen" */
bool
end_cond_sampling_Nielsen(const TwoStateCTMarkovModel &the_model,
                          const size_t start_state, const size_t end_state,
                          const double time_interval,
                          std::mt19937 &gen, std::vector<double> &jump_times,
                          const double start_time = 0.0);

/* "unif" */
void
end_cond_sample_unif(const TwoStateCTMarkovModel &the_model,
                     const size_t start_state, const size_t end_state,
                     const double time_interval,
                     std::mt19937 &gen, std::vector<double> &jump_times,
                     const double start_time = 0.0);

/* "pois" */
void
end_cond_sample_Poisson(const TwoStateCTMarkovModel &the_model,
                        const size_t start_state, const size_t end_state,
                        const double time_interval,
                        std::mt19937 &gen, std::vector<double> &jump_times,
                        const double start_time = 0.0);

size_t
forward_sampling(std::vector<std::function<double()> > &the_distrs,
                 size_t a, const double T, const double start_time,
                 std::vector<double> &jump_times);


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// #####  #####   ####  #####   ####   ###   #   ##     #####  #####   ####  #####
// ##  ## ##  ## ##  ## ##  ## ##  ## #     ###  ##     ##  ## ##  ## ##  ## ##  ##
// #####  #####  ##  ## #####  ##  ## #### #   # ##     #####  #####  ##  ## #####
// ##     ## ##  ##  ## ##     ##  ##    # ##### ##     ##     ## ##  ##  ## ##  ##
// ##     ##  ##  ####  ##      ####  ###  #   # #####  ##     ##  ##  ####  #####

double
end_cond_sample_prob(const TwoStateCTMarkovModel &the_model,
                     const std::vector<double> &jump_times,
                     const size_t start_state, const size_t end_state,
                     const double start_time, const double end_time,
                     size_t start_jump, const size_t end_jump);

#endif
