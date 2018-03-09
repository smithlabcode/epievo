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

void
model_param::get_branch_lengths(vector<double> &branch_lengths) const {
  t.get_branch_lengths(branch_lengths);
}


////////////////////////////////////////////////////////////////////////////////
// stationary and Markov chain
////////////////////////////////////////////////////////////////////////////////

void
potential_to_transprob(const vector<vector<double> > &Q,
                       vector<vector<double> > &T) {
  const double delta = sqrt(pow(Q[0][0] - Q[1][1], 2) + 4*Q[0][1]*Q[1][0]);
  // transition probability matrix
  T = Q;

  // compute the diagonal entries
  const double diag_denom = Q[0][0]+Q[1][1] + delta;
  T[1][1] = 2*Q[1][1]/diag_denom;
  T[0][0] = 2*Q[0][0]/diag_denom;

  // now compute the anti-diagonal entries
  const double anti_numer = 4*Q[0][1]*Q[1][0];
  T[0][1] = anti_numer/(pow(Q[0][0] + delta, 2) - Q[1][1]*Q[1][1]);
  T[1][0] = anti_numer/(pow(Q[1][1] + delta, 2) - Q[0][0]*Q[0][0]);
}



void
convert_parameter(const vector<vector<double> > &stationary_logfac,
                  vector<vector<double> > &T) {
  vector<vector<double> > Q = stationary_logfac;
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      Q[i][j] = exp(stationary_logfac[i][j]);

  potential_to_transprob(Q, T);
}


double
rate_factor(const vector<double> &rates) {
  // pair-wise potentials Q
  vector<vector<double> > Q(2, vector<double>(2, 1.0));
  Q[0][0] = Q[0][1] * sqrt(rates[2]/rates[0]);
  Q[1][1] = Q[0][1] * sqrt(rates[2]/rates[0]) * (rates[1]/rates[3]);

  // transition probability matrix
  vector<vector<double> > T;
  potential_to_transprob(Q, T);

  vector<double> pi(2, 0);
  pi[1] = (1.0 - T[0][0])/(2.0 - T[0][0] - T[1][1]);
  pi[0] = 1.0 - pi[1];
  vector<double> stationary_prob(rates.size(), 0.0);
  double unit = 0.0;
  for (size_t i = 0; i < rates.size(); ++i) {
    const size_t left = i/4;
    const size_t mid = (i % 4)/2;
    const size_t right = i % 2;
    stationary_prob[i] = pi[left]*T[left][mid]*T[mid][right];
    unit += stationary_prob[i] * rates[i];
  }

  return unit;
}
