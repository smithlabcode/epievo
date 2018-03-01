#ifndef JUMP_HPP
#define JUMP_HPP

#include <string>
#include <vector>

#include "path.hpp"
#include "TripletPattern.hpp"

struct Jump {
  double t;
  size_t pos;
  size_t context;
  std::vector<size_t> freq;
  Jump(double time, size_t position): t(time), pos(position) {};
};

void
get_jumps(const std::vector<Path> &paths, std::vector<Jump> &jumps);

void
get_scaled_jumps(const std::vector<Path> &paths, std::vector<Jump> &jumps);

void
get_suff_stat(const std::vector<std::vector<Jump> > &jumps,
              std::vector<double> &J, std::vector<double> &D);

void
get_suff_stat_by_branch(const std::vector<std::vector<Jump> > &jumps,
                        std::vector<double> &J,
                        std::vector<std::vector<double> > &D);

#endif
