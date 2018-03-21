/* Copyright (C) 2018 University of Southern California
 *                    Jianghan Qu and Andrew D Smith
 *
 * Author: Andrew D. Smith and Jianghan Qu
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA
 */

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <bitset>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "PhyloTreePreorder.hpp"
#include "Path.hpp"  /* related to Path */

#include "EpiEvoModel.hpp" /* model_param */
#include "TripletPattern.hpp"
// #include "jump.hpp"
#include "StateSeq.hpp"

#include "SingleSampler.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::min;
using std::runtime_error;
using std::bitset;

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


static void
propose_path(vector<std::exponential_distribution<double> > &distrs,
             std::mt19937 &gen,
             const EpiEvoModel &the_model,
             const Path &left, const Path &right, Path &mid) {

  size_t triplet = triple2idx(left.init_state, mid.init_state, right.init_state);
  double current_time = 0.0;

  auto i = left.jumps.begin();
  auto k = right.jumps.begin();

  mid.jumps.clear();

  // assume left and right have their final jump set to the total
  // time; we can fix this later
  while (i != left.jumps.end() && k != right.jumps.end()) {
    const double next_context_change = min(*i, *k);

    // sample a holding time in the current state
    const double holding = distrs[triplet](gen);

    if (current_time + holding < next_context_change) {
      current_time += holding;             // update the timepoint
      mid.jumps.push_back(current_time);   // add the timepoint to mid jumps
      triplet = flip_mid_bit(triplet);     // update the context
    }
    else {
      current_time = next_context_change;
      if (*i < *k) {
        triplet = flip_left_bit(triplet);
        ++i;
      }
      else {
        triplet = flip_right_bit(triplet);
        ++k;
      }
    }
    cerr << endl;
  }
}

int main(int argc, const char **argv) {
  try {

    bool VERBOSE = false;

    size_t the_site = 0;
    string node_name;

    ////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]), "test triple path",
                           "<model-file> <paths-file> <outfile>");
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("site", 's', "site to simulate",
                      true, the_site);
    opt_parse.add_opt("node", 'n', "name of node below edge to sample",
                      true, node_name);
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
    if (leftover_args.size() != 3) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string model_file(leftover_args[0]);
    const string pathsfile(leftover_args[1]);
    const string outfile(leftover_args[2]);
    ////////////////////////////////////////////////////////////////////////

    if (VERBOSE)
      cerr << "[READING PATHS: " << pathsfile << "]" << endl;
    vector<vector<Path> > all_paths; // along multiple branches
    vector<string> node_names;
    read_paths(pathsfile, node_names, all_paths);

    const size_t n_nodes = node_names.size();
    /* below: 1st element of all_paths empty at root; use last */
    const size_t n_sites = all_paths.back().size();

    if (VERBOSE)
      cerr << "n_nodes=" << n_nodes << endl
           << "n_sites=" << n_sites << endl;

    if (VERBOSE)
      cerr << "[READING PARAMETER FILE: " << model_file << "]" << endl;
    EpiEvoModel the_model;
    read_model(model_file, the_model);
    if (VERBOSE)
      cerr << the_model << endl;

    /* standard mersenne_twister_engine seeded with rd()*/
    std::random_device rd;
    std::mt19937 gen(rd());

    // get the id for the desired node
    vector<string>::const_iterator name_idx =
      find(begin(node_names), end(node_names), node_name);
    if (name_idx == node_names.end())
      throw runtime_error("invalid node name: " + node_name);

    const size_t node_id = name_idx - node_names.begin();
    const size_t parent_id = the_model.parent_ids[node_id];
    const double branch_length = the_model.branches[node_id];

    if (VERBOSE) {
      cerr << "node name: " << node_name << endl
           << "node id: " << node_id << endl
           << "parent name: " << node_names[parent_id] << endl
           << "parent id: " << parent_id << endl
           << "branch length: " << branch_length << endl
           << "site: " << the_site << endl
           << "total jumps: " << all_paths[node_id][the_site].jumps.size() << endl;
    }

    vector<std::exponential_distribution<double> > distrs;
    for (size_t i = 0; i < the_model.triplet_rates.size(); ++i) {
      auto ed = std::exponential_distribution<double>(the_model.triplet_rates[i]);
      cout << bitset<3>(i) << '\t' << ed << endl;
      distrs.push_back(ed);
    }

    Path left_path(all_paths[node_id][the_site - 1]);
    left_path.jumps.push_back(left_path.tot_time);

    Path right_path(all_paths[node_id][the_site + 1]);
    right_path.jumps.push_back(right_path.tot_time);

    Path the_path(all_paths[node_id][the_site]);
    propose_path(distrs, gen, the_model, left_path, right_path, the_path);

    cout << all_paths[node_id][the_site].end_state() << endl;
    cout << the_path.end_state() << endl;



    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    cerr << "------------ TEST homogeneous interval end-conditioned sampling BELOW -------------" << endl;
    double T = 0.5;
    vector<vector<double> > Q(2, vector<double>(2, 0.0));
    Q[0][1] = the_model.triplet_rates[0]; // 000
    Q[0][0] = - Q[0][1];
    Q[1][0] = the_model.triplet_rates[2]; // 010
    Q[1][1] = - Q[1][0];
    size_t a = 1;
    size_t b = 1;
    cerr << "Setup: interval of time\t " << T << endl
         << "       neighbor context\t 0-1" << endl
         << "       start state     \t " << a << endl
         << "       end state       \t " << b << endl;


    vector<double> jump_times;
    end_cond_sample(Q, a,  b, T, gen, jump_times);

    cerr << " sampled jump times: " << endl;
    for (size_t i = 0; i < jump_times.size(); ++i)
      cerr << jump_times[i] << ",\t" ;
    cerr << endl << "total " << jump_times.size()-2 << " jumps" << endl;;



  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
