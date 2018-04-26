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
  void get_branch_lengths(std::vector<double> &branch_lengths) const;
};


void
potential_to_transprob(const std::vector<std::vector<double> > &Q,
                       std::vector<std::vector<double> > &T);

void
convert_parameter(const std::vector<std::vector<double> > &stationary_logfac,
                  std::vector<std::vector<double> > &T);

double
rate_factor(const std::vector<double> &rates);

#endif
