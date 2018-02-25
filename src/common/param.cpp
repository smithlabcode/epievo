#include "param.hpp"

#include <fstream>
#include <sstream>

#include <cassert>
#include <cmath>

using std::vector;
using std::string;

void
model_param::read_param(const string &param_file) {
  std::ifstream in(param_file.c_str());
  if (!in)
    throw std::runtime_error("Could not open file" + param_file);

  string dummy_label;
  in >> dummy_label >> t;
  in >> dummy_label >> n_site;

  stationary_logfac =
    vector<vector<double> >(2, vector<double>(2, 0.0));
  stationary_logbaseline =
    vector<vector<double> >(2, vector<double>(2, 0.0));
  init_logfac =
    vector<vector<double> >(2, vector<double>(2, 0.0));

  in >> dummy_label >> stationary_logfac[0][0]
     >> stationary_logfac[0][1] >> stationary_logfac[1][1];
  assert(dummy_label == "stationary");
  stationary_logfac[1][0] = stationary_logfac[0][1];

  in >> dummy_label >> stationary_logbaseline[0][0]
     >> stationary_logbaseline[0][1] >> stationary_logbaseline[1][1];
  assert(dummy_label == "baseline");
  stationary_logbaseline[1][0] = stationary_logbaseline[0][1];

  in >> dummy_label >> init_logfac[0][0]
     >> init_logfac[0][1] >> init_logfac[1][1];
  assert(dummy_label == "init");
  init_logfac[1][0] = init_logfac[0][1];
}

void
model_param::get_rates(vector<double> &rates) const {
  rates.resize(8, 0);
  for (size_t i = 0; i < 2; ++i) {
    for (size_t j = 0; j < 2; ++j) {
      for (size_t k = 0; k < 2; ++k) {
        rates[4*i+2*j+k] = std::exp(stationary_logfac[i][1-j] +
                                    stationary_logfac[1-j][k] +
                                    stationary_logbaseline[i][k]);
      }
    }
  }
}
