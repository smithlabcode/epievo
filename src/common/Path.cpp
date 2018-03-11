#include "path.hpp"

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cassert>
#include <algorithm>    // std::lower_bound, std::upper_bound, std::sort

#include "smithlab_utils.hpp"

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

void initialize_paths(const std::vector<bool> &seq, const double tot_time,
                      std::vector<Path> &paths) {
  paths.resize(seq.size());
  for (size_t i = 0; i < seq.size(); ++i) {
    paths[i].init_state = seq[i];
    paths[i].tot_time = tot_time;
    paths[i].jumps.resize(0);
  }
}

void end_sequence(const vector<Path> &paths,
                  vector<bool> &seq) {
  seq.resize(paths.size());
  for (size_t i = 0; i < paths.size(); ++i) {
    bool s = paths[i].init_state;
    seq[i] = (paths[i].jumps.size() % 2 == 0)? s : !s;
  }
}

bool state_at_time(const Path &p, const double t) {
  vector<double>::const_iterator low;
  low = std::lower_bound(p.jumps.begin(), p.jumps.end(), t);
  size_t idx = (size_t)(low- p.jumps.begin());
  bool s = (idx%2 == 0)?  p.init_state: !p.init_state;
  // bool s = p.init_state;
  // for (size_t i = 0; i < p.jumps.size() && p.jumps[i] < t; ++i)
  //   s = !s;
  return s;
}

void sequence_at_time(const vector<Path> &paths, const double t,
                      vector<bool> &seq) {
  seq.resize(paths.size());
  for (size_t i = 0; i < paths.size(); ++i) {
    seq[i] = state_at_time(paths[i], t);
  }
}

void to_path(const bool s, const string jumps, Path &p) {
  p.init_state = s;
  vector<string> fields = smithlab::split(jumps, ",", false);
  p.tot_time = std::stod(fields.back());
  p.jumps.clear();
  for (size_t i = 1; i < fields.size()-1; ++i)
    p.jumps.push_back(std::stod(fields[i]));
}

void read_paths(const string path_file, vector<vector<Path> > &paths) {
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

void get_initial_seq(const vector<Path> &paths, vector<bool> &seq) {
  seq.clear();
  for (size_t i = 0; i < paths.size(); ++i)
    seq.push_back(paths[i].init_state);
}

////////////////////////////////////////////////////////////////////////////////

void intersect_paths(const Path &pa, const Path &pb, Environment &env) {
  assert(pa.tot_time == pb.tot_time);
  bool sa = pa.init_state;
  bool sb = pb.init_state;
  size_t i = 0;
  size_t j = 0;
  env.tot_time = pa.tot_time;
  while (i < pa.jumps.size() || j < pb.jumps.size()) {
    env.left.push_back(sa);
    env.right.push_back(sb);
    if (i < pa.jumps.size() && j < pb.jumps.size()) {
      if (pa.jumps[i] < pb.jumps[j]) {
        env.breaks.push_back(pa.jumps[i]);
        ++i;
        sa = !sa;
      } else if (pa.jumps[i] > pb.jumps[j]) {
        env.breaks.push_back(pb.jumps[j]);
        ++j;
        sb = !sb;
      } else {
        env.breaks.push_back(pb.jumps[j]);
        ++j;
        ++i;
        sa = !sa;
        sb = !sb;
      }
    } else if (i < pa.jumps.size()) {
      env.breaks.push_back(pa.jumps[i]);
      ++i;
      sa = !sa;
    } else {
      env.breaks.push_back(pb.jumps[j]);
      ++j;
      sb = !sb;
    }
  }

  if (env.breaks.size() == 0 ||
      env.breaks.back() < env.tot_time) {
    env.left.push_back(sa);
    env.right.push_back(sb);
  }
}

////////////////////////////////////////////////////////////////////////////////

TriplePath::TriplePath(const Path &l, const Path &m, const Path &r) {
  assert(l.tot_time == m.tot_time && l.tot_time == r.tot_time);
  size_t sl = (size_t)(l.init_state);
  size_t sm = (size_t)(m.init_state);
  size_t sr = (size_t)(r.init_state);
  states.clear();
  states.push_back(sl * 4 + sm * 2 + sr);

  breaks.clear();
  breaks.insert(breaks.end(), l.jumps.begin(), l.jumps.end());
  breaks.insert(breaks.end(), m.jumps.begin(), m.jumps.end());
  breaks.insert(breaks.end(), r.jumps.begin(), r.jumps.end());
  std::sort(breaks.begin(), breaks.end());
  breaks.push_back(l.tot_time);

  for (size_t i = 1; i < breaks.size(); ++i) {
    const double t = breaks[i-1] + (breaks[i] - breaks[i-1])/2;
    const size_t state_l = (size_t)(state_at_time(l, t));
    const size_t state_m = (size_t)(state_at_time(m, t));
    const size_t state_r = (size_t)(state_at_time(r, t));
    states.push_back(state_l * 4 + state_m * 2 + state_r);
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
