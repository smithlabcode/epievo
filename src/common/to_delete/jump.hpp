#ifndef JUMP_HPP
#define JUMP_HPP

#include <string>
#include <vector>

#include "Path.hpp"
#include "TripletPattern.hpp"

struct Jump {
  double t;
  size_t pos;
  size_t context;
  std::vector<size_t> freq;
  Jump(const double time, const size_t position): t(time), pos(position) {};
};

struct Hold {
  double hold_time;
  std::vector<size_t> freq;
};

void
get_jumps(const std::vector<Path> &paths, std::vector<Jump> &jumps,
          Hold &hold);

void
get_suff_stat(const std::vector<std::vector<Jump> > &jumps,
              const std::vector<Hold> &holds,
              std::vector<double> &J, std::vector<double> &D);

#endif
