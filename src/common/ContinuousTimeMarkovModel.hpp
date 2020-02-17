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
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA
 */

#ifndef CONTINUOUS_TIME_MARKOV_MODEL_HPP
#define CONTINUOUS_TIME_MARKOV_MODEL_HPP

#include <vector>
#include <string>

typedef std::vector<std::vector<double> > two_by_two;

struct CTMarkovModel {
  CTMarkovModel(const std::vector<double> &rates);
  CTMarkovModel(const std::pair<double, double> &rates);
  CTMarkovModel(const double r0, const double r1);
  void get_trans_prob_mat(const double time_interval,
                          two_by_two &prob_mat) const;
  std::string tostring() const;

  double get_rate(const bool the_state) const {
    return the_state ? rate1 : rate0;
  }
  double get_rate(const size_t the_state) const {
    return (the_state == 0ul) ? rate0 : rate1;
  }

  double rate0;
  double rate1;
  two_by_two U;
  two_by_two Uinv;
  std::vector<double> eigen_values;
};

std::ostream &
operator<<(std::ostream &os, const CTMarkovModel &ctmm);

struct TwoStateCTMarkovModel {
  TwoStateCTMarkovModel(const double r0, const double r1) :
    rate0(r0), rate1(r1) {}
  double get_trans_prob(const double time_interval, const size_t start_state,
                        const size_t end_state) const;
  std::string tostring() const;

  double get_rate(const bool the_state) const {
    return the_state ? rate1 : rate0;
  }
  double get_rate(const size_t the_state) const {
    return (the_state == 0ul) ? rate0 : rate1;
  }

  double rate0;
  double rate1;
};

std::ostream &
operator<<(std::ostream &os, const TwoStateCTMarkovModel &ctmm);


void
continuous_time_trans_prob_mat(const double rate0, const double rate1,
                               const double time_interval,
                               two_by_two &transition_matrix);

////////////////////////////////////////////////////////////////////////////////
///////////////////       CTMM-RELATED STATISTICS     //////////////////////////
////////////////////////////////////////////////////////////////////////////////

void
expectation_J(const double r0, const double r1, const double T,
              two_by_two &J0, two_by_two &J1);

void
expectation_D(const double r0, const double r1, const double T,
              two_by_two &D0, two_by_two &D1);

#endif
