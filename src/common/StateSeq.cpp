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
using std::runtime_error;

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
                 triplet_props.begin(),
                 std::bind(std::divides<double>(), std::placeholders::_1,
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
  std::transform(pair_counts.begin(), pair_counts.end(), pair_props.begin(),
                 std::bind(std::divides<double>(), std::placeholders::_1,
                           seq.size() - 1));
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

std::ostream &
operator<<(std::ostream &os, const StateSeq &s) {
  for (size_t i = 0; i < s.seq.size(); ++i)
    os << (s.seq[i] == true ? '1' : '0');
  return os;
}

#include <sstream>
#include <fstream>

void
read_states_file(const string &statesfile,
                 vector<string> &node_names,
                 vector<StateSeq> &the_states) {

  std::ifstream in(statesfile.c_str());
  if (!in)
    throw runtime_error("cannot read states file: " + statesfile);

  string buffer;
  getline(in, buffer);
  if (buffer[0] == '#') // ADS: probably should adopt some convention
    buffer = buffer.substr(1);

  std::istringstream iss(buffer);
  string tmp_name;
  while (iss >> tmp_name)
    node_names.push_back(tmp_name);

  const size_t n_nodes = node_names.size();

  size_t the_site = 0; // dummy
  the_states.resize(n_nodes);
  while (getline(in, buffer)) {
    iss.clear();
    iss.str(std::move(buffer));

    iss >> the_site;

    size_t val_count = 0;
    char tmp_value = 0;
    while (val_count < n_nodes && iss >> tmp_value)
      the_states[val_count++].seq.push_back(tmp_value == '1');

    if (val_count != n_nodes)
      throw runtime_error("bad line in states file");
  }
}
