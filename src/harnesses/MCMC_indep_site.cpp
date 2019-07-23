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

using std::function;
using std::exponential_distribution;
using std::placeholders::_1;


////////////////////////////////////////////////////////////////////////////////

template <class T>
static void
print_vec(vector<T> vec) {
  for (size_t i = 0; i < vec.size(); i++)
    cerr << vec[i] << ", ";
  cerr << endl;
}

static void
collect_init_sequences(const vector<Path> paths, vector<bool> &sequences) {
  sequences.clear();
  for (size_t site_id = 0; site_id < paths.size(); site_id++)
    sequences.push_back(paths[site_id].init_state == 1);
}

static void
collect_end_sequences(const vector<Path> paths, vector<bool> &sequences) {
  sequences.clear();
  for (size_t site_id = 0; site_id < paths.size(); site_id++)
    sequences.push_back(paths[site_id].end_state() == 1);
}

static void
get_indep_params(const EpiEvoModel &the_model, vector<double> &rates,
                 double &pi_0) {
  
  rates.resize(2, 0.0);
  
  vector<double> tri_rates = the_model.triplet_rates;
  
  rates[0] = (tri_rates[0] + tri_rates[1] + tri_rates[4] + tri_rates[5]) / 4;
  rates[1] = (tri_rates[2] + tri_rates[3] + tri_rates[6] + tri_rates[7]) / 4;
  
  pi_0 = the_model.init_T[0][0] / (the_model.init_T[0][0] + the_model.init_T[1][1]);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
static void
sample_root(const double pi_0, const size_t n_sites,
            vector<bool> &seq, std::mt19937 &gen) {

  seq.resize(n_sites, false);
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  
  for (size_t i = 0; i < n_sites; i++)
    seq[i] = (unif(gen) > pi_0);
}



static size_t
forward_sampling(vector<function<double()> > &the_distrs,
                 size_t a, const double T,
                 vector<double> &jump_times) {
  jump_times.clear();
  double tau = 0.0;
  while ((tau += the_distrs[a]()) < T) {
    a = complement_state(a);
    jump_times.push_back(tau);
  }
  return a;
}


static bool
check_leaf_seq(const vector<bool> seq, const vector<Path> &by_site) {

  bool pass = true;
  size_t site_id = 0;
  while(pass && site_id < by_site.size()) {
    pass = (by_site[site_id].end_state() == seq[site_id]) ? true : false;
    site_id++;
  }
  return pass;
}


static void
forward_simulation(const TwoStateCTMarkovModel &ctmm, const double pi_0,
                   const TreeHelper &th,
                   const size_t n_sites, vector<vector<Path> > &paths,
                   std::mt19937 &gen) {
  
  const size_t n_nodes = th.n_nodes;
  paths.resize(n_nodes);

  vector<bool> root_seq;
  sample_root(pi_0, n_sites, root_seq, gen);

  typedef exponential_distribution<double> exp_distr;
  vector<function<double()> > the_distrs = {
    function<double()>(bind(exp_distr(ctmm.rate0), ref(gen))),
    function<double()>(bind(exp_distr(ctmm.rate1), ref(gen)))
  };
  
  for (size_t site_id = 0; site_id < n_sites; site_id++)
    paths[0].push_back(Path(root_seq[site_id], 0.0));
  
  for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
    vector<Path> path_by_site(n_sites);
    const double curr_branch_len = th.branches[node_id];

    for (size_t site_id = 0; site_id < n_sites; site_id++) {
      const size_t start_state = paths[th.parent_ids[node_id]][site_id].end_state();
      path_by_site[site_id].init_state = start_state;
      path_by_site[site_id].tot_time = curr_branch_len;

      forward_sampling(the_distrs, start_state, curr_branch_len,
                       path_by_site[site_id].jumps);
    }
    paths[node_id] = path_by_site;
  }
}


static bool
end_cond_forward_simulation(const TwoStateCTMarkovModel &ctmm,
                            const double pi_0, const TreeHelper &th,
                            const vector<vector<bool>> end_sequences,
                            vector<vector<Path> > &paths, std::mt19937 &gen,
                            const bool FIX_ROOT) {
  
  const size_t n_nodes = end_sequences.size();
  const size_t n_sites = end_sequences[0].size();
  
  vector<bool> root_seq(end_sequences[0]);
  if (!FIX_ROOT)
    sample_root(pi_0, n_sites, root_seq, gen);

  for (size_t site_id = 0; site_id < n_sites; site_id++)
    paths[0].push_back(Path(root_seq[site_id], 0.0));
  
  typedef exponential_distribution<double> exp_distr;
  vector<function<double()> > the_distrs = {
    function<double()>(bind(exp_distr(ctmm.rate0), ref(gen))),
    function<double()>(bind(exp_distr(ctmm.rate1), ref(gen)))
  };

  bool pass = true;
  for (size_t node_id = 1; pass && (node_id < n_nodes); ++node_id) {
    // start new branch
    const double curr_branch_len = th.branches[node_id];
    vector<Path> path_by_site(n_sites);
    
    for (size_t site_id = 0; site_id < n_sites; site_id++) {
      const size_t start_state = paths[th.parent_ids[node_id]][site_id].end_state();
      path_by_site[site_id].init_state = start_state;
      path_by_site[site_id].tot_time = curr_branch_len;
      const size_t end_state = end_sequences[node_id][site_id];

      forward_sampling(the_distrs, start_state, curr_branch_len,
                       path_by_site[site_id].jumps);
    }

    paths[node_id] = path_by_site;
    
    // leaf sequence agree?
    pass = th.is_leaf(node_id) ?
    check_leaf_seq(end_sequences[node_id], paths[node_id]) : true;
  }
  return pass;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* generate initial paths by heuristics */
static void
initialize_paths(std::mt19937 &gen, const TreeHelper &th,
                 vector<vector<bool> > &end_sequences,
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
///////////////////               PRUNING              /////////////////////////
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
double
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
bool
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



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

static void
write_statistics_header(const string outfile) {
  std::ofstream out(outfile.c_str());
  out << "ITR\tSAMPLE\tJ_0\tJ_1\tD_0\tD_1" << endl;
}


static void
write_statistics(const string outfile, const size_t itr, const size_t sample,
                 const vector<double> &J, const vector<double> &D) {
  std::ofstream out(outfile.c_str(), std::ofstream::app);
  out << itr << "\t" << sample << "\t";
  for (size_t i = 0; i < 2; i++)
    out << J[i] << "\t";
  for (size_t i = 0; i < 1; i++)
    out << D[i] << "\t";
  out << D.back() << endl;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char **argv) {
  try {
    bool VERBOSE = false;
    //bool FIX_ROOT = false;
    
    string outfile;
    string tree_file;
    
    size_t n_sites = 5;
    size_t batch = 100;
    size_t n_mcmc_batches = 100;

    size_t rng_seed = std::numeric_limits<size_t>::max();
    
    double evolutionary_time = numeric_limits<double>::lowest();
    
    ///////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]),
                           "test MCMC against forward-simulation summary stats",
                           "site-independent <param>");
    opt_parse.add_opt("n_sites", 'n', "number of sites", false, n_sites);
    opt_parse.add_opt("mcmc_itr", 'i', "number of MCMC iterations",
                      false, n_mcmc_batches);
    opt_parse.add_opt("batch", 'B', "batch size",
                      false, batch);
    opt_parse.add_opt("tree", 't', "Newick format tree file", false, tree_file);
    opt_parse.add_opt("evo-time", 'T', "evolutionary time", false,
                      evolutionary_time);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);;
    opt_parse.add_opt("outfile", 'o', "outfile (prefix)",
                      false, outfile);
    //opt_parse.add_opt("fixroot", 'R', "fix root states",
    //                  false, FIX_ROOT);

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
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    const string param_file(leftover_args[0]);
    ///////////////////////////////////////////////////////////////////////////

    /* (1) READING PARAMETERS FROM FILE */
    if (VERBOSE)
      cerr << "[READING PARAMETERS: " << param_file << "]" << endl;
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
    
    if (VERBOSE)
      cerr << the_model << endl;

    /* (2) LOADING TREE */
    size_t n_nodes = 0;
    TreeHelper th;
    if (!tree_file.empty()) {
      if (VERBOSE)
        cerr << "[READING TREE: " << tree_file << "]" << endl;
      PhyloTreePreorder the_tree; // tree topology and branch lengths
      std::ifstream tree_in(tree_file.c_str());
      if (!tree_in || !(tree_in >> the_tree))
        throw std::runtime_error("bad tree file: " + tree_file);
      n_nodes = the_tree.get_size();
      th = TreeHelper(the_tree);
    }
    else {
      if (VERBOSE)
        cerr << "[INITIALIZING TWO-NODE TREE WITH TIME: "
        << evolutionary_time << "]" << endl;
      n_nodes = 2;
      th = TreeHelper(evolutionary_time);
    }
     
    /* (3) INITIALIZING THE RANDOM NUMBER GENERATOR */
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    std::mt19937 gen(rng_seed);
    if (VERBOSE)
      cerr << "[INITIALIZED RNG (seed=" << rng_seed << ")]" << endl;

    /* ********************************************************************** */
    /* (4) GENERATE ONE PATH */
    if (VERBOSE)
      cerr << "[GENERATE A SAMPLE PATH]" << endl;
    vector<vector<Path> > paths; // along multiple branches
    forward_simulation(ctmm, pi_0, th, n_sites, paths, gen);
  
    // collect sequences
    vector<bool> root_seq;
    collect_init_sequences(paths[1], root_seq);
    vector<vector<bool>> end_sequences (n_nodes, root_seq);
    for (size_t node_id = 1; node_id < n_nodes; node_id++) {
      collect_end_sequences(paths[node_id], end_sequences[node_id]);
    }
    /* ********************************************************************** */
    /* (5) FORWARD SIMULATION */
    if (VERBOSE)
      cerr << "[FORWARD SIMULATION]" << endl;

    // forward-simulation output files
    vector<string> statfiles_forward(n_nodes);
    for (size_t node_id = 1; node_id < n_nodes; node_id++) {
      const string fstat = outfile + "." + th.node_names[node_id] + ".forward";
      statfiles_forward[node_id] = fstat;
      write_statistics_header(fstat);
    }

    // sampling
    size_t n_forward_samples_collected = 0;
    while (n_forward_samples_collected < batch) {
      vector<vector<Path> > sampled_paths (n_nodes);

      const bool success = end_cond_forward_simulation(ctmm, pi_0, th,
                                                       end_sequences,
                                                       sampled_paths, gen,
                                                       true);
      if (success) {
        vector<vector<double> > J, D;
        compute_sufficient_statistics(sampled_paths, J, D);
        for (size_t node_id = 1; node_id < n_nodes; node_id++)
          write_statistics(statfiles_forward[node_id],
                           0, n_forward_samples_collected,
                           J[node_id], D[node_id]);
        n_forward_samples_collected++;
      }
    }

    /* ********************************************************************** */
    /* (6) MCMC */

    // mcmc output files
    vector<string> statfiles_mcmc(n_nodes);
    for (size_t node_id = 1; node_id < n_nodes; node_id++) {
      const string fstat = outfile + "." + th.node_names[node_id] + ".mcmc";
      statfiles_mcmc[node_id] = fstat;
      write_statistics_header(fstat);
    }
    
    // distort/randomize paths
    vector<vector<Path> > mcmc_paths(paths);
    initialize_paths(gen, th, end_sequences, mcmc_paths, true);

    vector<Path> proposed_path;
    for(size_t batch_id = 0; batch_id < n_mcmc_batches; batch_id++) {
      for (size_t sample_id = 0; sample_id < batch; sample_id++) {
        for (size_t site_id = 0; site_id < n_sites; ++site_id) {
          MCMC_indep_site(ctmm, pi_0, th, site_id, mcmc_paths, gen,
                          proposed_path);
          assert(mcmc_paths[1][site_id].end_state() == paths[1][site_id].end_state());
        }
        
        // write stats
        vector<vector<double> > J, D;
        compute_sufficient_statistics(mcmc_paths, J, D);
        for (size_t node_id = 1; node_id < n_nodes; node_id++)
          write_statistics(statfiles_mcmc[node_id],
                           batch_id, sample_id, J[node_id], D[node_id]);
      }
    }
    
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
