#ifndef PARAM_HPP
#define PARAM_HPP

#include <string>
#include <vector>

#include "PhyloTreePreorder.hpp"

struct model_param {
  PhyloTreePreorder t;
  size_t n_site;
  std::vector<std::vector<double> > stationary_logfac; // 2 by 2
  std::vector<std::vector<double> > stationary_logbaseline;  // 2 by 2
  std::vector<std::vector<double> > init_logfac;

  void read_param(const std::string &param_file);
  void get_rates(std::vector<double> &rates) const;
};


#endif
