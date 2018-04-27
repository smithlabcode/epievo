/* Copyright (C) 2018 University of Southern California
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

#include "Path.hpp"

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cassert>
#include <algorithm>    // std::lower_bound, std::upper_bound, std::sort

#include "smithlab_utils.hpp"

#include "StateSeq.hpp"

using std::vector;
using std::string;

static const string NODE_TAG = "NODE";
static const size_t TAG_LENGTH = 4;

static bool
is_node_line(const string &buffer) {
  return buffer.length() > TAG_LENGTH &&
    (buffer.compare(0, TAG_LENGTH, NODE_TAG) == 0);
}

static string
get_node_name(const string &buffer) {
  return string(buffer.begin() + buffer.find(':') + 1, buffer.end());
}

std::ostream &
operator<<(std::ostream &os, const Path &p) {
  os << p.init_state << '\t' << p.tot_time << '\t';
  copy(p.jumps.begin(), p.jumps.end(), std::ostream_iterator<double>(os, "\t"));
  return os;
}

void
initialize_paths(const std::vector<bool> &seq, const double tot_time,
                 std::vector<Path> &paths) {
  paths.resize(seq.size());
  for (size_t i = 0; i < seq.size(); ++i)
    paths[i] = Path(seq[i], tot_time);
}

void
get_seq_end(const vector<Path> &paths, vector<bool> &seq) {
  seq.resize(paths.size());
  for (size_t i = 0; i < paths.size(); ++i)
    seq[i] = paths[i].end_state();
}

bool
Path::state_at_time(const double t) const {
  const size_t idx =
    std::lower_bound(jumps.begin(), jumps.end(), t) - jumps.begin();
  const bool s = (idx % 2 == 0) ? init_state : !init_state;
  return s;
}

void
get_seq_at_time(const double t, const vector<Path> &paths,
                vector<bool> &seq) {
  seq.resize(paths.size());
  for (size_t i = 0; i < paths.size(); ++i) {
    seq[i] = paths[i].state_at_time(t);
  }
}

void
to_path(const bool s, const string &jumps, Path &p) {
  p.init_state = s;
  vector<string> fields = smithlab::split(jumps, ",", false);
  p.tot_time = std::stod(fields.back());
  p.jumps.clear();
  for (size_t i = 1; i < fields.size()-1; ++i)
    p.jumps.push_back(std::stod(fields[i]));
}

void
read_paths(const string path_file, vector<vector<Path> > &paths) {
  std::ifstream in(path_file.c_str());
  if (!in)
    throw SMITHLABException("cannot read: " + path_file);

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
    throw std::runtime_error("cannot read: " + path_file);

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

void get_seq_init(const vector<Path> &paths, vector<bool> &seq) {
  seq.clear();
  for (size_t i = 0; i < paths.size(); ++i)
    seq.push_back(paths[i].init_state);
}

////////////////////////////////////////////////////////////////////////////////

Environment::Environment(const Path &pa, const Path &pb) {
  assert(pa.tot_time == pb.tot_time);
  bool sa = pa.init_state;
  bool sb = pb.init_state;
  size_t i = 0;
  size_t j = 0;
  tot_time = pa.tot_time;
  while (i < pa.jumps.size() || j < pb.jumps.size()) {
    left.push_back(sa);
    right.push_back(sb);
    if (i < pa.jumps.size() && j < pb.jumps.size()) {
      if (pa.jumps[i] < pb.jumps[j]) {
        breaks.push_back(pa.jumps[i]);
        ++i;
        sa = !sa;
      } else if (pa.jumps[i] > pb.jumps[j]) {
        breaks.push_back(pb.jumps[j]);
        ++j;
        sb = !sb;
      } else {
        breaks.push_back(pb.jumps[j]);
        ++j;
        ++i;
        sa = !sa;
        sb = !sb;
      }
    } else if (i < pa.jumps.size()) {
      breaks.push_back(pa.jumps[i]);
      ++i;
      sa = !sa;
    } else {
      breaks.push_back(pb.jumps[j]);
      ++j;
      sb = !sb;
    }
  }

  if (breaks.size() == 0 ||
      breaks.back() < tot_time) {
    left.push_back(sa);
    right.push_back(sb);
    breaks.push_back(tot_time);
  }
}

////////////////////////////////////////////////////////////////////////////////

TriplePath::TriplePath(const Path &l, const Path &m, const Path &r) {
  assert(l.tot_time == m.tot_time && l.tot_time == r.tot_time);
  states.clear();
  states.push_back(triple2idx(l.init_state, m.init_state, r.init_state));

  breaks.clear();
  breaks.insert(breaks.end(), l.jumps.begin(), l.jumps.end());
  breaks.insert(breaks.end(), m.jumps.begin(), m.jumps.end());
  breaks.insert(breaks.end(), r.jumps.begin(), r.jumps.end());
  // consider two uses of std::inplace_merge below, rather than sort
  std::sort(breaks.begin(), breaks.end());
  breaks.push_back(l.tot_time);

  for (size_t i = 1; i < breaks.size(); ++i) {
    const double t = breaks[i-1] + (breaks[i] - breaks[i-1])/2;
    states.push_back(triple2idx(l.state_at_time(t),
                                m.state_at_time(t), r.state_at_time(t)));
  }

  jump_context_freq.resize(8, 0);
  for (size_t i = 0; i < m.jumps.size(); ++i) {
    vector<double>::iterator low;
    low = std::lower_bound (breaks.begin(), breaks.end(), m.jumps[i]);
    const size_t context = states[(size_t)(low - breaks.begin())];
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


PathContextStat::PathContextStat(const Path &l, const Path &m, const Path &r) {
  jumps_in_context = vector<double>(8, 0.0);
  time_in_context = vector<double>(8, 0.0);
  size_t context = triple2idx(l.init_state, m.init_state, r.init_state);
  vector<double> jumps = m.jumps;
  jumps.insert(jumps.begin(), 0.0);
  jumps.push_back(m.tot_time);
  for (size_t i = 1; i < jumps.size(); ++i) {
    const double t = 0.5*(jumps[i] + jumps[i-1]);
    context = triple2idx(l.state_at_time(t), m.state_at_time(t),
                         r.state_at_time(t));
    ++jumps_in_context[context];
    time_in_context[context] += jumps[i] - jumps[i-1];
  }
  --jumps_in_context[context]; // last break point is not a jump
}
