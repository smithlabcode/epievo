#ifndef PATH_HPP
#define PATH_HPP

#include <string>
#include <vector>
#include <fstream>
#include <sstream>

using std::vector;
using std::string;

struct Path {
  bool init_state;
  double tot_time;
  vector<double> jumps;
};

struct Environment {
  vector<bool> left;
  vector<bool> right;
  vector<double> breaks;
  double tot_time;
};

void initialize_paths(const std::vector<bool> &seq, const double tot_time,
                      std::vector<Path> &paths) {
  paths.resize(seq.size());
  for (size_t i = 0; i < seq.size(); ++i) {
    paths[i].init_state = seq[i];
    paths[i].tot_time = tot_time;
    paths[i].jumps.resize(0);
  }
}

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

void end_sequence(const vector<Path> &paths,
                  vector<bool> &seq) {
  seq.resize(paths.size());
  for (size_t i = 0; i < paths.size(); ++i) {
    bool s = paths[i].init_state;
    seq[i] = (paths[i].jumps.size() % 2 == 0)? s : !s;
  }
}

bool state_at_time(const Path &p, const double t) {
  bool s = p.init_state;
  for (size_t i = 0; i < p.jumps.size() && p.jumps[i] < t; ++i)
    s = !s;
  return s;
}

void sequence_at_time(const vector<Path> &paths,
                      const double t,
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

void get_inital_seq(const vector<Path> paths, vector<bool> &seq) {
  seq.clear();
  for (size_t i = 0; i < paths.size(); ++i)
    seq.push_back(paths[i].init_state);
}


#endif
