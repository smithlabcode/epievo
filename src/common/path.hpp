#ifndef PATH_HPP
#define PATH_HPP

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cassert>

#include "smithlab_utils.hpp"

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
                      std::vector<Path> &paths);

void intersect_paths(const Path &pa, const Path &pb, Environment &env);

void end_sequence(const vector<Path> &paths, vector<bool> &seq);

bool state_at_time(const Path &p, const double t);

void sequence_at_time(const vector<Path> &paths,
                      const double t, vector<bool> &seq);

void to_path(const bool s, const string jumps, Path &p);

void read_paths(const string path_file, vector<vector<Path> > &paths);

void get_inital_seq(const vector<Path> paths, vector<bool> &seq);

#endif
