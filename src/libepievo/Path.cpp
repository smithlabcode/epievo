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

#include "Path.hpp"

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cassert>
#include <iomanip>
#include <algorithm>
#include <exception>
#include <iterator>

#include "epievo_utils.hpp"

using std::vector;
using std::string;
using std::to_string;
using std::runtime_error;
using std::ostream_iterator;
using std::istringstream;

static const string NODE_TAG = "NODE";
static const size_t TAG_LENGTH = 4;

static bool
is_node_line(const string &buffer) {
  return buffer.length() > TAG_LENGTH &&
    (buffer.compare(0, TAG_LENGTH, NODE_TAG) == 0);
}

static string
get_node_name(const string &buffer) {
  return string(begin(buffer) + buffer.find(':') + 1, end(buffer));
}

std::ostream &
operator<<(std::ostream &os, const Path &p) {
  std::ios old_state(nullptr);
  old_state.copyfmt(os);
  os << p.init_state << '\t' << p.tot_time << '\t';
  copy(begin(p.jumps), end(p.jumps), ostream_iterator<double>(os, "\t"));
  os.precision(std::numeric_limits<double>::max_digits10);
  os.copyfmt(old_state);
  return os;
}

void
initialize_paths(const state_seq &seq,
                 const double tot_time, vector<Path> &paths) {
  paths.resize(seq.size());
  for (size_t i = 0; i < seq.size(); ++i)
    paths[i] = Path(seq[i], tot_time);
}

void
get_seq_end(const vector<Path> &paths, state_seq &seq) {
  seq.resize(paths.size());
  for (size_t i = 0; i < paths.size(); ++i)
    seq[i] = paths[i].end_state();
}

void
Path::scale_to_unit_length() {
  transform(begin(jumps), end(jumps), begin(jumps),
            [&](const double d) {return d/tot_time;});
  tot_time = 1.0;
}

bool
Path::state_at_time(const double t) const {
  const size_t idx =
    std::lower_bound(begin(jumps), end(jumps), t) - begin(jumps);
  const bool s = (idx % 2 == 0) ? init_state : !init_state;
  return s;
}

void
get_seq_at_time(const double t, const vector<Path> &paths,
                state_seq &seq) {
  seq.resize(paths.size());
  for (size_t i = 0; i < paths.size(); ++i) {
    seq[i] = paths[i].state_at_time(t);
  }
}

void
to_path(const bool s, const string &jumps, Path &p) {
  p.init_state = s;
  istringstream iss(jumps);
  string item;
  getline(iss, item, ','); // eat the first one
  p.jumps.clear();
  while (getline(iss, item, ','))
    p.jumps.push_back(std::stod(item));
  p.tot_time = p.jumps.back();
  p.jumps.resize(p.jumps.size() - 1);
}

void
read_paths(const string path_file, vector<vector<Path> > &paths) {
  std::ifstream in(path_file.c_str());
  if (!in)
    throw runtime_error("cannot read: " + path_file);

  vector<size_t> node_ids;
  size_t node_id, pos;
  bool init_state;
  string jumpstring;
  string line;
  while (std::getline(in, line)) {
    if (line[0]!= '#') {
      std::istringstream iss(line);
      iss >> node_id >> pos >> init_state >> jumpstring;
      if (node_ids.size() == 0 || node_id != node_ids.back()) {
        node_ids.push_back(node_id);
        vector<Path> pp;
        paths.push_back(pp);
      }
      Path p;
      to_path(init_state, jumpstring, p);
      paths[node_ids.size() - 1].push_back(p);
    }
  }
}


void
read_paths(const string &path_file, vector<string> &node_names,
           vector<vector<Path> > &paths) {

  std::ifstream in(path_file.c_str());
  if (!in)
    throw runtime_error("cannot read: " + path_file);

  size_t pos = 0;
  double tmp_jump = 0.0;
  string line;
  while (getline(in, line)) {
    if (is_node_line(line)) {
      node_names.push_back(get_node_name(line));
      paths.push_back(vector<Path>());
    }
    else {
      Path p;
      std::istringstream iss(std::move(line));
      iss >> pos >> p.init_state >> p.tot_time;
      while (iss >> tmp_jump)
        p.jumps.push_back(tmp_jump);
      paths.back().push_back(p);
    }
  }
}

void get_seq_init(const vector<Path> &paths, state_seq &seq) {
  seq.clear();
  for (size_t i = 0; i < paths.size(); ++i)
    seq.push_back(paths[i].init_state);
}

////////////////////////////////////////////////////////////////////////////////

struct TriplePath {
  vector<size_t> states; // triplet states, length k
  vector<double> breaks; // start is first jump, end is total_time, length k
  vector<size_t> jump_context_freq; // context freq. of jumps at middle site

  TriplePath(const Path &l, const Path &m, const Path &r);
  void time_by_context(vector<double> &tbc) const;
};

TriplePath::TriplePath(const Path &l, const Path &m, const Path &r) {
  assert(l.tot_time == m.tot_time && l.tot_time == r.tot_time);
  states.clear();
  states.push_back(triple2idx(l.init_state, m.init_state, r.init_state));

  breaks.clear();
  breaks.insert(end(breaks), begin(l.jumps), end(l.jumps));
  breaks.insert(end(breaks), begin(m.jumps), end(m.jumps));
  breaks.insert(end(breaks), begin(r.jumps), end(r.jumps));
  // consider two uses of std::inplace_merge below, rather than sort
  std::sort(begin(breaks), end(breaks));
  breaks.push_back(l.tot_time);

  for (size_t i = 1; i < breaks.size(); ++i) {
    // ADS: this, below, is bad...
    const double t = breaks[i-1] + (breaks[i] - breaks[i-1])/2;
    states.push_back(triple2idx(l.state_at_time(t),
                                m.state_at_time(t),
                                r.state_at_time(t)));
  }

  jump_context_freq.resize(8, 0);
  for (size_t i = 0; i < m.jumps.size(); ++i) {
    vector<double>::iterator low;
    low = std::lower_bound (begin(breaks), end(breaks), m.jumps[i]);
    const size_t context = states[(size_t)(low - begin(breaks))];
    ++jump_context_freq[context];
  }
}

void
TriplePath::time_by_context(vector<double> &tbc) const {
  tbc.resize(8, 0.0);
  tbc[states[0]] += breaks[0];
  for (size_t i = 0; i < states.size(); ++i)
    tbc[states[i]] += breaks[i+1] - breaks[i];
}


static void
add_suff_stat_single(const size_t change_pos, const Path &a, size_t i,
                     size_t triplet, double prev_time,
                     vector<double> &J, vector<double> &D) {
  while (i < a.jumps.size()) {
    D[triplet] += a.jumps[i] - prev_time;
    // only update J if jumping in middle path
    if (change_pos == 1) J[triplet] += 1.0;
    triplet = ((change_pos == 1) ? flip_mid_bit(triplet) :
               (change_pos == 0) ? flip_left_bit(triplet) :
               flip_right_bit(triplet));
    prev_time = a.jumps[i];
    ++i;
  }
  // need to count the final interval towards D but never count J for
  // final segment
  D[triplet] += a.tot_time - prev_time;
}

/* Process a pair of paths collecting sufficient statistics J and
   D. The pair (a, b) may be a (left,mid), (left,righ) or (mid,right).
   The "fixed" position parameter indicates which member of the
   original three paths is not present, and therefore which bit
   doesn't change in the triplet. The original positions of a and b
   determines which bit is flipped in the triplet at the end of
   segments, and if J is updated (for a, for b, or neither). Once one
   member of the pair has no more jumps, the jumps in the other member
   of the pair are processed as a "single" using add_suff_stat_single,
   which has a parameter for the "changing" position, instead of the
   "fixed" position.
 */
static void
add_suff_stat_pair(const size_t fixed_pos,
                   const Path &a, size_t i, const Path &b, size_t j,
                   size_t triplet, double prev_time,
                   vector<double> &J, vector<double> &D) {
  assert(fixed_pos <= 2);

  while (i < a.jumps.size() && j < b.jumps.size())
    if (a.jumps[i] < b.jumps[j]) {
      D[triplet] += a.jumps[i] - prev_time;
      if (fixed_pos == 0) // no change at left, "a" is mid, count J
        J[triplet] += 1.0;
      triplet = fixed_pos == 0 ? flip_mid_bit(triplet) : flip_left_bit(triplet);
      prev_time = a.jumps[i];
      ++i;
    }
    else {
      D[triplet] += b.jumps[j] - prev_time;
      if (fixed_pos == 2) // no change at right, "b" is mid, count J
        J[triplet] += 1.0;
      triplet = fixed_pos == 2 ? flip_mid_bit(triplet) : flip_right_bit(triplet);
      prev_time = b.jumps[j];
      ++j;
    }
  if (i < a.jumps.size())
    add_suff_stat_single(fixed_pos == 0 ? 1 : 0, a, i, triplet, prev_time, J, D);
  else
    add_suff_stat_single(fixed_pos == 2 ? 1 : 2, b, j, triplet, prev_time, J, D);
}

void
add_sufficient_statistics(const Path &left, const Path &mid, const Path &right,
                          vector<double> &J, vector<double> &D) {

  size_t triplet = triple2idx(left.init_state, mid.init_state, right.init_state);
  double prev_time = 0.0;
  size_t i = 0, j = 0, k = 0;
  while (i < left.jumps.size() && j < mid.jumps.size() && k < right.jumps.size())
    if (left.jumps[i] < std::min(mid.jumps[j], right.jumps[k])) { // LEFT
      D[triplet] += left.jumps[i] - prev_time;
      prev_time = left.jumps[i];
      triplet = flip_left_bit(triplet);
      ++i;
    }
    else if (mid.jumps[j] < right.jumps[k]) { // MID
      D[triplet] += mid.jumps[j] - prev_time;
      J[triplet] += 1.0; // update J because we are at mid position
      prev_time = mid.jumps[j];
      triplet = flip_mid_bit(triplet);
      ++j;
    }
    else { // RIGHT (also: default case, shouldn't udpate J)
      D[triplet] += right.jumps[k] - prev_time;
      prev_time = right.jumps[k];
      triplet = flip_right_bit(triplet);
      ++k;
    }

  if (i == left.jumps.size())
    add_suff_stat_pair(0, mid, j, right, k, triplet, prev_time, J, D);
  else if (j == mid.jumps.size())
    add_suff_stat_pair(1, left, i, right, k, triplet, prev_time, J, D);
  else // (k == right.jumps.size())
    add_suff_stat_pair(2, left, i, mid, j, triplet, prev_time, J, D);
}
