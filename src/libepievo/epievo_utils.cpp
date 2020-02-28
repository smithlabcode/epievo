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

#include "epievo_utils.hpp"

#include <string>
#include <vector>
#include <sstream>
#include <functional>
#include <algorithm>
#include <fstream>
#include <stdexcept>

using std::string;
using std::vector;
using std::runtime_error;

void
get_triplet_counts(const state_seq &seq, vector<size_t> &triplet_counts) {
  static const size_t n_triplets = 8;
  triplet_counts.resize(n_triplets, 0);
  for (size_t i = 2; i < seq.size(); ++i)
    triplet_counts[triple2idx(seq[i-2], seq[i-1], seq[i])]++;
}

void
get_triplet_proportions(const state_seq &seq, vector<double> &triplet_props) {
  vector<size_t> triplet_counts;
  get_triplet_counts(seq, triplet_counts);
  triplet_props.resize(triplet_counts.size());
  const size_t denom = seq.size() - 2;
  transform(begin(triplet_counts), end(triplet_counts),
            begin(triplet_props), [denom](const double x) {return x/denom;});
}

void
get_pair_counts(const state_seq &seq, vector<size_t> &pair_counts) {
  static const size_t n_pairs = 4;
  pair_counts.resize(n_pairs, 0);
  for (size_t i = 1; i < seq.size(); ++i)
    pair_counts[pair2idx(seq[i-1], seq[i])]++;
}

void
get_pair_proportions(const state_seq &seq,
                     vector<double> &pair_props) {
  vector<size_t> pair_counts;
  get_pair_counts(seq, pair_counts);
  pair_props.resize(pair_counts.size());
  const size_t denom = seq.size() - 1;
  transform(begin(pair_counts), end(pair_counts), begin(pair_props),
            [denom](const double x) {return x/denom;});
}

string
summary_string(const state_seq &s) {

  vector<double> triplet_props;
  get_triplet_proportions(s, triplet_props);
  vector<double> pair_props;
  get_pair_proportions(s, pair_props);

  std::ostringstream oss;
  oss << "triplet proportions:\n"
      << triplet_info_to_string(triplet_props) << '\n'
      << "pair proportions:\n"
      << pair_info_to_string(pair_props);
  return oss.str();
}

void
read_states_file(const string &statesfile, vector<string> &names,
                 vector<state_seq> &the_states) {

  std::ifstream in(statesfile);
  if (!in)
    throw runtime_error("cannot read states file: " + statesfile);

  string buffer;
  getline(in, buffer);
  if (buffer[0] == '#') // ADS: probably should adopt some convention
    buffer = buffer.substr(1);

  std::istringstream iss(buffer);
  string tmp_name;
  while (iss >> tmp_name)
    names.push_back(tmp_name);

  const size_t n_seqs = names.size();

  size_t the_site = 0; // dummy
  the_states.resize(n_seqs);
  while (getline(in, buffer)) {
    iss.clear();
    iss.str(std::move(buffer));

    iss >> the_site;

    size_t val_count = 0;
    char tmp_value = 0;
    while (val_count < n_seqs && iss >> tmp_value)
      the_states[val_count++].push_back(tmp_value == '1');

    if (val_count != n_seqs)
      throw std::runtime_error("bad line in states file");
  }
}
