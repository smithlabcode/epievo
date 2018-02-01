#ifndef PARAM_HPP
#define PARAM_HPP

#include <string>
#include <vector>
#include <fstream>
#include <sstream>

using std::vector;
using std::string;

struct model_param {
  PhyloTreePreorder t;
  size_t n_site;
  // double tot_time;
  vector<vector<double> > stationary_logfac; // 2 by 2
  vector<vector<double> > stationary_logbaseline;  // 2 by 2
  vector<vector<double> > init_logfac;
};

void read_param(const string param_file, model_param &p) {
  std::ifstream in(param_file.c_str());
  if (!in)
    throw std::runtime_error("Could not open file" + param_file);

  string dummy_label;
  in >> dummy_label >> p.t;
  in >> dummy_label >> p.n_site;
  //  in >> dummy_label >> p.tot_time;
  p.stationary_logfac =
    vector<vector<double> >(2, vector<double>(2, 0.0));
  p.stationary_logbaseline =
    vector<vector<double> >(2, vector<double>(2, 0.0));
  p.init_logfac =
    vector<vector<double> >(2, vector<double>(2, 0.0));
  in >> dummy_label >> p.stationary_logfac[0][0]
     >> p.stationary_logfac[0][1] >> p.stationary_logfac[1][1];
  assert(dummy_label == "stationary");
  p.stationary_logfac[1][0] = p.stationary_logfac[0][1];
  in >> dummy_label >> p.stationary_logbaseline[0][0]
     >> p.stationary_logbaseline[0][1] >> p.stationary_logbaseline[1][1];
  assert(dummy_label == "baseline");
  p.stationary_logbaseline[1][0] = p.stationary_logbaseline[0][1];
  in >> dummy_label >> p.init_logfac[0][0]
     >> p.init_logfac[0][1] >> p.init_logfac[1][1];
  assert(dummy_label == "init");
  p.init_logfac[1][0] = p.init_logfac[0][1];
}

void get_rates(const model_param &p, vector<double> &rates) {
  rates.resize(8, 0);
  for (size_t i = 0; i < 2; ++i) {
    for (size_t j = 0; j < 2; ++j) {
      for (size_t k = 0; k < 2; ++k) {
       rates[4*i+2*j+k] = exp(p.stationary_logfac[i][1-j] +
                              p.stationary_logfac[1-j][k] +
                              p.stationary_logbaseline[i][k]);
      }
    }    
  }  
}


#endif
