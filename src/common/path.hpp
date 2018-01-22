#ifndef PATH_HPP
#define PATH_HPP

#include <string>
#include <vector>

using std::vector;

struct path {
  bool init_state;
  double tot_time;
  vector<double> jumps;
};

struct environment {
  vector<bool> left;
  vector<bool> right;
  vector<double> breaks;
  double tot_time;
};

void initialize_paths(const std::vector<bool> &seq, const double tot_time,
                        std::vector<path> &paths) {
  paths.resize(seq.size());
  for (size_t i = 0; i < seq.size(); ++i) {
    paths[i].init_state = seq[i];
    paths[i].tot_time = tot_time;
    paths[i].jumps.resize(0);
  }
}

void intersect_paths(const path &pa, const path &pb, environment &env) {
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

void end_sequence(const vector<path> &paths,
                  vector<bool> &seq) {
  seq.resize(paths.size());
  for (size_t i = 0; i < paths.size(); ++i) {
    bool s = paths[i].init_state;
    seq[i] = (paths[i].jumps.size() % 2 == 0)? s : !s;
  }
}

bool state_at_time(const path &p, const double t) {
  bool s = p.init_state;
  for (size_t i = 0; i < p.jumps.size() && p.jumps[i] < t; ++i)
    s = !s;
  return s;
}

void sequence_at_time(const vector<path> &paths,
                      const double t,
                      vector<bool> &seq) {
  seq.resize(paths.size());
  for (size_t i = 0; i < paths.size(); ++i) {
    seq[i] = state_at_time(paths[i], t);
  }
}

#endif
