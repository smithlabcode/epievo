/* Copyright (C) 2019 University of Southern California
 *                    Xiaojing Ji, Jianghan Qu and Andrew D Smith
 *
 * Author: Andrew D. Smith, Jianghan Qu and Xiaojing Ji
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
#include <sstream>
#include <algorithm>
#include <bitset>
#include <random>
#include <functional>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "Path.hpp"
#include "EpiEvoModel.hpp"
#include "StateSeq.hpp"
#include "EndCondSampling.hpp"
#include "ContinuousTimeMarkovModel.hpp"
#include "TreeHelper.hpp"
#include "IndepSite.hpp"
#include "TripletSampler.hpp"
#include "GlobalJump.hpp"
#include "ParamEstimation.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::runtime_error;
using std::numeric_limits;
using std::to_string;

///////////////////////////////////////////////////////////////////////////
static void
get_indep_params(const EpiEvoModel &the_model, vector<double> &rates,
                 double &pi_0) {
  
  rates.resize(2, 0.0);
  
  vector<double> tri_rates = the_model.triplet_rates;
  
  rates[0] = (tri_rates[0] + tri_rates[1] + tri_rates[4] + tri_rates[5]) / 4;
  rates[1] = (tri_rates[2] + tri_rates[3] + tri_rates[6] + tri_rates[7]) / 4;
  
  pi_0 = the_model.init_T[0][0] / (the_model.init_T[0][0] + the_model.init_T[1][1]);
}


template <class T>
size_t
read_states_file(const string &statesfile, vector<vector<T> > &state_sequences,
                 TreeHelper &th) {
  
  std::ifstream in(statesfile.c_str());
  if (!in)
    throw std::runtime_error("bad states file: " + statesfile);
  
  // first line is the list of node names
  string buffer;
  if (!getline(in, buffer))
    throw std::runtime_error("cannot read nodes line in: " + statesfile);
  
  std::istringstream nodes_iss(buffer);
  vector<string> node_names;
  
  string tmp_node_name;
  while (nodes_iss >> tmp_node_name)
    node_names.push_back(tmp_node_name);
  
  if (node_names.size() != 2)
    throw std::runtime_error("more or fewer than 2 nodes in: " + statesfile);
  
  if (node_names.front()[0] == '#') {
    if (node_names.front().length() == 1)
      node_names = vector<string>(node_names.begin() + 1, node_names.end());
    else node_names.front() = node_names.front().substr(1);
  }
  
  th.node_names = node_names;
  
  // read states
  const size_t n_nodes = node_names.size();
  size_t site_count = 0;
  
  state_sequences.resize(th.n_nodes);
  
  while (getline(in, buffer)) {
    std::istringstream iss;
    iss.str(std::move(buffer));
    
    size_t site_index = 0;
    iss >> site_index; // not important info but must be removed
    
    // now read the states for the current site
    size_t node_idx = 0;
    T tmp_state_val;
    while (node_idx < n_nodes && iss >> tmp_state_val)
        state_sequences[node_idx++].push_back(tmp_state_val);
    
    if (node_idx < n_nodes)
      throw std::runtime_error("inconsistent number of states: " +
                               to_string(node_idx) + "/" + to_string(n_nodes));
    ++site_count;
  }
  
  if (site_count == 0)
    throw std::runtime_error("no sites read from states file: " + statesfile);

  return site_count;
}



/* generate initial paths by heuristics */
static void
initialize_paths_indep(std::mt19937 &gen, const TreeHelper &th,
                       vector<vector<size_t> > &end_sequences,
                       vector<vector<Path> > &paths, const bool FIX_ROOT) {
  
  const size_t n_sites = end_sequences.front().size();
  
  paths.resize(th.n_nodes);
  
  auto unif =
  bind(std::uniform_real_distribution<double>(0.0, 1.0), std::ref(gen));
  
  vector<size_t> child_states = {0, 0}; // two child states
  
  for (size_t i = th.n_nodes; i > 0; --i) {
    
    const size_t node_id = i - 1;
    paths[node_id].resize(n_sites);
    
    for (size_t site_id = 0; site_id < n_sites; ++site_id) {
      
      if (!th.is_leaf(node_id)) {
        // get child states
        size_t n_ch = 0;
        for (ChildSet c(th.subtree_sizes, node_id); c.good(); ++c)
          child_states[n_ch++] = end_sequences[*c][site_id];
        
        // sample parent state
        if (th.parent_ids[node_id] > 0 || !FIX_ROOT)
          end_sequences[node_id][site_id] = child_states[std::floor(unif()*n_ch)];
        
        // assign site-specific paths above two child nodes
        for (ChildSet c(th.subtree_sizes, node_id); c.good(); ++c) {
          
          const double len = th.branches[*c];
          paths[*c][site_id].tot_time = len;
          paths[*c][site_id].init_state = end_sequences[node_id][site_id];
          
          const bool parent_state = end_sequences[node_id][site_id];
          if (parent_state != paths[*c][site_id].init_state)
            paths[*c][site_id].jumps.push_back(unif() * len);
        }
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////                 MCMC               /////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* Iterates over all nodes, in reverse pre-order, and for each it sets
 * the "p" and "q" from Felsenstein's algorithm, requiring for each
 * node that the children have already been processed (hence reversing
 * the order).
 */
static void
pruning_upward(const TreeHelper &th, const size_t site_id,
               const vector<vector<Path> > &paths,
               const TwoStateCTMarkovModel &ctmm, vector<FelsHelper> &fh) {
  
  fh.resize(th.n_nodes);
  for (size_t i = th.n_nodes; i > 0; --i) {
    const size_t node_id = i - 1;
    
    vector<double> q(2, 1.0);
    if (th.is_leaf(node_id)) {
      const bool leaf_state = paths[node_id][site_id].end_state();
      q[0] = (leaf_state == false) ? 1.0 : 0.0;
      q[1] = (leaf_state == true)  ? 1.0 : 0.0;
    }
    else {
      for (ChildSet c(th.subtree_sizes, node_id); c.good(); ++c) {
        q[0] *= fh[*c].p[0];
        q[1] *= fh[*c].p[1];
      }
    }
    fh[node_id].q.swap(q); // assign computed q to the fh
    
    /* now calculate p if we are not at root */
    if (!th.is_root(node_id)) {
      vector<vector<double> > P;
      continuous_time_trans_prob_mat(ctmm.rate0, ctmm.rate1,
                                     paths[node_id][site_id].tot_time, P);
      // p <- P*q
      vector<double> p = { P[0][0]*fh[node_id].q[0] + P[0][1]*fh[node_id].q[1],
        P[1][0]*fh[node_id].q[0] + P[1][1]*fh[node_id].q[1] };
      fh[node_id].p.swap(p); // assign computed p to the fh
    }
  }
}


/* Iterates over all nodes in pre-order, and for each branch it samples
 * a new path.
 */
static void
sampling_downward(const TwoStateCTMarkovModel &ctmm, const double &root_p0,
                  const TreeHelper &th, const size_t site_id,
                  const vector<vector<Path> > &paths,
                  std::mt19937 &gen, vector<Path> &proposed_path,
                  double &orig_proposal, double &update_proposal,
                  const vector<FelsHelper> &fh, const bool FIX_ROOT) {
  
  // compute posterior probability at root node
  proposed_path = vector<Path>(th.n_nodes);
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  proposed_path.front() = Path(paths[1][site_id].init_state, 0.0); // copy root
  if (!FIX_ROOT)
    proposed_path.front() = Path(unif(gen) > root_p0, 0.0); // update root
  
  update_proposal = proposed_path.front().init_state ?
  log(1.0 - root_p0) : log(root_p0);
  orig_proposal = paths[1][site_id].init_state ?
  log(1.0 - root_p0) : log(root_p0);
  
  for (size_t node_id = 1; node_id < th.n_nodes; ++node_id) {
    const size_t start_state = proposed_path[th.parent_ids[node_id]].end_state();
    const double T = paths[node_id][site_id].tot_time;
    
    proposed_path[node_id].init_state = start_state;
    proposed_path[node_id].tot_time = T;
    
    vector<vector<double> > P;
    continuous_time_trans_prob_mat(ctmm.rate0, ctmm.rate1, T, P);
    const double p0 =
    P[start_state][0] * fh[node_id].q[0] / fh[node_id].p[start_state];
    const size_t sampled_state = (unif(gen) > p0);
    end_cond_sample_Poisson(ctmm, start_state, sampled_state,
                            T, gen, proposed_path[node_id].jumps, 0.0);
    orig_proposal +=
    end_cond_sample_Poisson_prob(ctmm, paths[node_id][site_id].jumps,
                                 start_state, paths[node_id][site_id].end_state(),
                                 0, T, 0, paths[node_id][site_id].jumps.size());
    update_proposal +=
    end_cond_sample_Poisson_prob(ctmm, proposed_path[node_id].jumps,
                                 start_state, sampled_state,
                                 0, T,
                                 0, proposed_path[node_id].jumps.size());
    assert(paths[node_id][site_id].end_state()==proposed_path[node_id].end_state());
  }
}


static double
indep_log_likelihood(const vector<double> &rates,
                     const vector<double> &J, const vector<double> &D) {
  double r = 0.0;
  for (size_t i = 0; i < rates.size(); ++i) {
    r += J[i]*log(rates[i]) - D[i]*rates[i];
  }
  return r;
}



static void
get_suff_stats(const vector<Path> &paths, vector<double> &J, vector<double> &D) {
  
  J.resize(2, 0.0);
  D.resize(2, 0.0);
  
  for (size_t b = 1; b < paths.size(); b++) {
    
    size_t prev_state = paths[b].init_state;
    double time = 0.0;
    
    for (size_t j = 0; j < paths[b].jumps.size(); j++) {
      J[prev_state] += 1;
      D[prev_state] += ( paths[b].jumps[j] - time );
      prev_state = 1 - prev_state;
      time = paths[b].jumps[j];
    }
    D[prev_state] += (paths[b].tot_time - time);
  }
}




/* compute acceptance rate */
static double
log_accept_rate(const vector<double> rates, const size_t site_id,
                const vector<vector<Path> > &paths,
                const vector<Path> &proposed_path,
                const double orig_proposal, const double update_proposal) {
  
  // calculate likelihood
  vector<Path> original(paths.size());
  for (size_t i = 1; i < paths.size(); ++i)
    original[i] = paths[i][site_id];
  
  double llr = orig_proposal - update_proposal;
  
  vector<double> J_orig, D_orig, J_prop, D_prop;
  get_suff_stats(original, J_orig, D_orig);
  get_suff_stats(proposed_path, J_prop, D_prop);
  
  llr += (indep_log_likelihood(rates, J_prop, D_prop) -
          indep_log_likelihood(rates, J_orig, D_orig));
  
  return llr;
}

////////////////////////////////////////////////////////////////////////////////
///////////////          single MCMC iteration         /////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* Metropolis-Hastings sampling at single site */
static bool
MCMC_indep_site(const TwoStateCTMarkovModel &ctmm, const double &pi_0,
                const TreeHelper &th,
                const size_t site_id, vector<vector<Path> > &paths,
                std::mt19937 &gen, vector<Path> &proposed_path) {
  
  vector<FelsHelper> fh;
  double orig_proposal, update_proposal;
  
  pruning_upward(th, site_id, paths, ctmm, fh);
  sampling_downward(ctmm, pi_0, th, site_id, paths, gen, proposed_path,
                    orig_proposal, update_proposal, fh, true);
  
  // acceptance rate
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  const double u = unif(gen) ;
  const double log_acc_rate =
  log_accept_rate({ctmm.rate0, ctmm.rate1}, site_id, paths, proposed_path,
                  orig_proposal, update_proposal);
  
  bool accepted = false;
  if (log_acc_rate >= 0 || u < exp(log_acc_rate))
    accepted = true;
  
  // if accepted, replace old path with proposed one.
  if (accepted)
    for (size_t i = 1; i < th.n_nodes; ++i)
      paths[i][site_id] = proposed_path[i];
  
  return accepted;
}


///////////////////////////////////////////////////////////////////////////

static void
write_root_to_pathfile_local(const string &outfile, const string &root_name) {
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream outpath(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  if (!outpath)
    throw std::runtime_error("bad output file: " + outfile);
  
  outpath << "NODE:" << root_name << endl;
}


static void
append_to_pathfile_local(const string &pathfile, const string &node_name,
                         const vector<Path> &path_by_site) {
  std::ofstream of;
  if (!pathfile.empty()) of.open(pathfile.c_str(), std::ofstream::app);
  std::ostream outpath(pathfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  if (!outpath)
    throw std::runtime_error("bad output file: " + pathfile);
  
  outpath << "NODE:" << node_name << endl;
  for (size_t i = 0; i < path_by_site.size(); ++i)
    outpath << i << '\t' << path_by_site[i] << '\n';
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char **argv) {
  try {
    bool VERBOSE = false;

    string outfile;               // sampled local path
    string pathfile;               // initial local path
    
    size_t burnin = 10;           // burn-in MCMC iterations
    double evolutionary_time = 1.0;

    size_t rng_seed = std::numeric_limits<size_t>::max();

    ///////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]),
                           "Simulate paths between a pair of sequences",
                           "<param> <states>");
    opt_parse.add_opt("burnin", 'L',
                      "MCMC burn-in length (default: 10)",
                      false, burnin);
    opt_parse.add_opt("evo-time", 'T', "evolutionary time", false,
                      evolutionary_time);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);
    opt_parse.add_opt("outfile", 'o',
                      "output file of local paths (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("pathfile", 'p',
                      "input file of local paths", false, pathfile);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);

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
    if (leftover_args.size() < 2) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    const string param_file(leftover_args[0]);
    const string states_file(leftover_args[1]);
    ///////////////////////////////////////////////////////////////////////////
    /* (1) READING PARAMETERS FROM FILE */
    if (VERBOSE)
      cerr << "[READING PARAMETERS: " << param_file << endl;
    EpiEvoModel the_model;
    read_model(param_file, the_model);
    the_model.scale_triplet_rates();
    vector<double> rates;
    double pi_0;
    get_indep_params(the_model, rates, pi_0);
    const TwoStateCTMarkovModel ctmm(rates[0], rates[1]);
    if (VERBOSE)
      cerr << "[R0: " << rates[0] << ", R1: " << rates[1]
      << ", PI0: " << pi_0 << "]" << endl;

    /* (2) GET THE ROOT SEQUENCE */
    if (VERBOSE)
      cerr << "[OBTAINING ROOT/LEAF SEQUENCE]" << endl;
    TreeHelper th(evolutionary_time);
    vector<vector<size_t> > state_sequences; // double if data is continuous
    read_states_file(states_file, state_sequences, th);
    const size_t n_sites = state_sequences.front().size();
    
    /* (3) INITIALIZING THE RANDOM NUMBER GENERATOR */
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    std::mt19937 gen(rng_seed);

    vector<vector<Path> > init_paths;
    vector<string> node_names;
    if (!pathfile.empty())
      read_paths(pathfile, node_names, init_paths);
    else 
      initialize_paths_indep(gen, th, state_sequences, init_paths, true);


    if (VERBOSE)
      cerr << "[SIMULATING]" << endl;
    
    
    /* INITIALIZING THE PATH */
    vector<vector<Path> > paths = init_paths;
    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    /////// STARTING MCMC
    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    
    /* METROPOLIS-HASTINGS ALGORITHM */
    // Burning
    double n_acc = 0;
    for(size_t burnin_itr = 0; burnin_itr < burnin; burnin_itr++) {
      for (size_t site_id = 1; site_id < n_sites - 1; ++site_id) {
        vector<Path> proposed_path;
        n_acc += MCMC_indep_site(ctmm, pi_0, th, site_id, paths, gen,
                                 proposed_path);
        assert(paths[1][site_id].init_state == state_sequences[0][site_id]);
        assert(paths[1][site_id].end_state() == state_sequences[1][site_id]);
      }
    }

    if (VERBOSE)
      cout << ", acc rate: " << n_acc / ((n_sites-2)*burnin) << endl;
    write_root_to_pathfile_local(outfile, th.node_names.front());
    for (size_t node_id = 1; node_id < th.n_nodes; ++node_id)
      append_to_pathfile_local(outfile, th.node_names[node_id], paths[node_id]);
  }

  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
