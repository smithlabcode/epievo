#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>    // std::min


#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "PhyloTreePreorder.hpp"
#include "Path.hpp"  /* related to Path */
#include "param.hpp" /* model_param */
#include "TripletPattern.hpp"
#include "jump.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::min;

void
propose_new_path(const vector<double> &rates, const Path &l,
                 const Path &m, const Path &r,
                 Path &new_m, std::mt19937 &gen) {
  size_t curr_context = l.init_state*4 + m.init_state*2 + r.init_state;
  double curr_t = 0;
  bool curr_state = m.init_state;
  new_m.init_state = m.init_state;
  new_m.tot_time = m.tot_time;
  new_m.jumps.clear();
  Environment env(l, r);
  for (size_t i = 1; i < env.breaks.size(); ++i) {
    curr_context = env.left[i]*4 + curr_state*2 + env.right[i];
    while (curr_t < env.breaks[i]) {
      std::exponential_distribution<double> exp_distr(rates[curr_context]);
      const double holding_time = exp_distr(gen);
      if (curr_t + holding_time < env.breaks[i]) {
        curr_t += holding_time;
        curr_state = !curr_state;
        curr_context = curr_context ^ 2ul; // flip middle bit
        new_m.jumps.push_back(curr_t);
      } else {
        curr_t = env.breaks[i];
      }
    }
  }
}

double
log_proposal_prob(const vector<double> &rates, const Path &l,
                  const Path &m, const Path &r) {
  size_t curr_context = l.init_state*4 + m.init_state*2 + r.init_state;
  double curr_t = 0;
  bool curr_state = m.init_state;
  Environment env(l, r);
  double log_prob = 0;
  size_t j = 0;
  for (size_t i = 1; i < env.breaks.size(); ++i) {
    while (j < m.jumps.size() && m.jumps[j] < env.breaks[i]) {
      curr_context = env.left[i]*4 + curr_state*2 + env.right[i];
      const double holding_time = m.jumps[j] - curr_t;
      log_prob += log(rates[curr_context]) - rates[curr_context]*holding_time;
      curr_t += holding_time;
      ++j;
      curr_state = !curr_state;
      curr_context = curr_context ^ 2ul; // flip middle bit
    }
    log_prob += -rates[curr_context]*(env.breaks[i] - curr_t);
    curr_t = env.breaks[i];
  }
  return log_prob;
}

static double
log_lik_ratio(const vector<double> &rates,
              const PathContextStat &pcs_num,
              const PathContextStat &pcs_denom) {
  double result = 0.0;
  for (size_t i = 0; i < 8; ++i) {
    result += (pcs_num.jumps_in_context[i] -
               pcs_denom.jumps_in_context[i])*log(rates[i]) -
      (pcs_num.time_in_context[i] -
       pcs_denom.time_in_context[i]) * rates[i];
  }
  return result;
}


double
log_accept_rate(const vector<double> &rates,
                const Path &ll, const Path &l, const Path &m,
                const Path &r, const Path &rr, const Path &new_m) {
  const double lp_prop_old = log_proposal_prob(rates, l, m, r);
  const double lp_prop_new = log_proposal_prob(rates, l, new_m, r);

  PathContextStat pcs_old(l, m, r);
  PathContextStat pcs_new(l, new_m, r);

  PathContextStat pcs_old_l(ll, l, m);
  PathContextStat pcs_new_l(ll, l, new_m);

  PathContextStat pcs_old_r(m, r, rr);
  PathContextStat pcs_new_r(new_m, r, rr);

  // cerr << "likratio_middle =" << log_lik_ratio(rates, pcs_new, pcs_old) << endl;
  // cerr << "likratio_left =" << log_lik_ratio(rates, pcs_new_l, pcs_old_l) << endl;
  // cerr << "likratio_right =" << log_lik_ratio(rates, pcs_new_r, pcs_old_r) << endl;

  double lr = lp_prop_old - lp_prop_new +
    log_lik_ratio(rates, pcs_new, pcs_old) +
    log_lik_ratio(rates, pcs_new_l, pcs_old_l) +
    log_lik_ratio(rates, pcs_new_r, pcs_old_r);
  return lr;
}

int main(int argc, const char **argv) {
  try {
    string param_file;
    string out_file;
    OptionParser opt_parse(strip_path(argv[0]), "test triple path",
                           "<path-file>");
    opt_parse.add_opt("param", 'p', "params file", true, param_file);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() < 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string path_file(leftover_args.back());

    cerr << "Reading parameters" << endl;
    model_param p;
    p.read_param(param_file);
    vector<double> rates;
    p.get_rates(rates);

    cerr << "Reading paths " << path_file << endl;
    vector<vector<Path> > paths; // along multiple branches
    read_paths(path_file, paths);

    ////////////////////////////////////////////////////////////
    // - paths to jumps
    // - collect suff_stats
    //
    vector<vector<Jump> > jumps;
    vector<Hold> holds;

    for (size_t b = 0; b < paths.size(); ++b) {
      vector<Jump> j;
      Hold h;
      get_jumps(paths[b], j, h);
      jumps.push_back(j);
      holds.push_back(h);
    }

    ////////////////////////////////////////////////////////////
    cerr << "Propose new path" << endl;
    /* standard mersenne_twister_engine seeded with rd()*/
    std::random_device rd;
    std::mt19937 gen(rd());

    for (size_t i = 0; i < paths.size(); ++i) {
      size_t tot_acc = 0;
      size_t diff = 0;
      for (size_t k = 1; k < paths[0].size()-3; ++k) {
        Path new_path;
        propose_new_path(rates, paths[0][k], paths[0][k+1], paths[0][k+2],
                         new_path, gen);

        double lar = log_accept_rate(rates, paths[0][k-1], paths[0][k],
                                     paths[0][k+1], paths[0][k+2],
                                     paths[0][k+3], new_path);

        std::uniform_real_distribution<double> unif(0.0,1.0);
        if (lar != 0) ++diff;
        if (lar != 0 && unif(gen) < min(1.0, exp(lar))) { //accept
          if (lar!= 0) ++tot_acc;
          paths[0][k+1] = new_path;
        }
      }
      cerr <<"branch-" << i << ":\tProposed "
           << diff << " new paths\tAccepted " << tot_acc << " paths \t"
           << "out of total "  << paths[0].size()-5 << " paths"<< endl;
    }
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
