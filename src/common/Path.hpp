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

std::ostream &
operator<<(std::ostream &os, const Path &p);

void to_path(const bool s, const string jumps, Path &p);

bool state_at_time(const Path &p, const double t);

void initialize_paths(const std::vector<bool> &seq, const double tot_time,
                      std::vector<Path> &paths);

void read_paths(const string path_file, vector<vector<Path> > &paths);

void
read_paths(const string &path_file,
           std::vector<std::string> &node_names,
           std::vector<std::vector<Path> > &paths);

void get_initial_seq(const vector<Path> &paths, vector<bool> &seq);

void end_sequence(const vector<Path> &paths, vector<bool> &seq);

void sequence_at_time(const vector<Path> &paths, const double t,
                      vector<bool> &seq);

////////////////////////////////////////////////////////////////////////////////

struct Environment {
  vector<bool> left; // states on the left
  vector<bool> right;  // states on the right
  vector<double> breaks; // time points of environment state changes
  double tot_time;
};


void intersect_paths(const Path &pa, const Path &pb, Environment &env);

////////////////////////////////////////////////////////////////////////////////
struct TriplePath {
  vector<size_t> states; /*triplet states, length k*/
  vector<double> breaks; /*start with first jump, end with total_time, length k*/
  vector<size_t> jump_context_freq; /*context frequency of jumps at the middle site*/

  TriplePath(const Path &l, const Path &m, const Path &r);
  void time_by_context(vector<double> &tbc) const;
};

#endif
