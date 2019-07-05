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
#include "SingleSiteSampler.hpp"
#include "IntervalSampler.hpp"
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


////////////////////////////////////////////////////////////////////////////////

static void
collect_init_sequences(const vector<Path> paths,
                                  StateSeq &sequences) {
  vector<char> seq;
  for (size_t site_id = 0; site_id < paths.size(); site_id++)
    seq.push_back('0' + paths[site_id].init_state);
  sequences.seq = seq;
}

static void
collect_end_sequences(const vector<Path> paths,
                                  StateSeq &sequences) {
  vector<char> seq;
  for (size_t site_id = 0; site_id < paths.size(); site_id++)
    seq.push_back('0' + paths[site_id].end_state());
  sequences.seq = seq;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* This function does the sampling for an individual change in the
state sequence
 */
static void
sample_jump(const EpiEvoModel &the_model, const double total_time,
            std::mt19937 &gen, TripletSampler &ts, vector<GlobalJump> &the_path,
            double &time_value) {
  static const size_t n_triplets = 8;
  
  // triplet_count = c_{ijk} for current sequence (encoded in the
  // TripletSampler object)
  vector<size_t> triplet_counts;
  ts.get_triplet_counts(triplet_counts);
  
  // holding_rate = c_{ijk}*lambda_{ijk}
  const double holding_rate =
  std::inner_product(triplet_counts.begin(), triplet_counts.end(),
                     the_model.triplet_rates.begin(), 0.0);
  
  // sample a holding time = time until next state change
  std::exponential_distribution<double> exp_distr(holding_rate);
  const double holding_time = std::max(exp_distr(gen),
                                       numeric_limits<double>::min());
  
  // update the current time_value
  time_value += holding_time;
  
  // if the holding time ends before the total time interval, we can
  // make a change to the state sequence
  if (time_value < total_time) {
    
    /* first: get a probability distribution for the triplet to change */
    vector<double> triplet_prob(n_triplets, 0.0);
    for (size_t i = 0; i < n_triplets; ++i)
      triplet_prob[i] =
      triplet_counts[i]*the_model.triplet_rates[i]/holding_rate;
    
    /* next: use that distribution to sample which triplet type at
     which the change will happen */
    std::discrete_distribution<size_t> multinom(triplet_prob.begin(),
                                                triplet_prob.end());
    const size_t context = multinom(gen);
    
    /* sample a change position having the relevant triplet; this
     changes the TripletSampler data structure to reflect a changed
     state at the position sampled */
    const size_t change_position = ts.random_mutate(context, gen);
    
    /* add the changed position and change time to the path */
    the_path.push_back(GlobalJump(time_value, change_position));
  }
}


static bool
check_leaf_seq(const StateSeq s, const vector<Path> &by_site) {
  bool pass = true;
  size_t site_id = 0;
  while(pass && site_id < by_site.size()) {
    const bool leaf_state = (s.seq[site_id] == '1');
    pass = (by_site[site_id].end_state() == leaf_state) ? true : false;
    site_id++;
  }
  return pass;
}


static void
assign_changes_to_sites(const vector<GlobalJump> &global_path,
                        vector<Path> &by_site) {
  
  const size_t n_changes = global_path.size();
  for (size_t i = 0; i < n_changes; ++i) {
    const size_t pos = global_path[i].position;
    const double time = global_path[i].timepoint;
    by_site[pos].jumps.push_back(time);
  }
}


static void
forward_simulation(const EpiEvoModel &the_model, const TreeHelper &th,
                   const size_t n_sites, vector<vector<Path> > &paths,
                   std::mt19937 &gen) {
  
  const size_t n_nodes = th.n_nodes;
  paths.resize(n_nodes);

  vector<char> root_seq;
  the_model.sample_state_sequence_init(n_sites, gen, root_seq);

  StateSeq s(root_seq);
  vector<StateSeq> sequences(n_nodes, s);
  
  for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
    // start new branch
    vector<GlobalJump> global_path;
    const double curr_branch_len = th.branches[node_id];
    
    TripletSampler ts(sequences[th.parent_ids[node_id]]);
    double time_value = 0;
    while (time_value < curr_branch_len)
      sample_jump(the_model, curr_branch_len, gen, ts, global_path,
                  time_value);
    ts.get_sequence(sequences[node_id]);
    
    // Convert global jumps to local paths
    vector<Path> path_by_site(n_sites);
    for (size_t site_id = 0; site_id < n_sites; site_id++) {
      path_by_site[site_id].init_state = s.seq[site_id];
      path_by_site[site_id].tot_time = curr_branch_len;
    }
    assign_changes_to_sites(global_path, path_by_site);
    paths[node_id] = path_by_site;
  }
}


static bool
end_cond_forward_simulation(const EpiEvoModel &the_model, const TreeHelper &th,
                            const vector<StateSeq> end_sequences,
                            vector<vector<Path> > &paths, std::mt19937 &gen,
                            const bool FIX_ROOT) {

  const size_t n_nodes = end_sequences.size();
  const size_t n_sites = end_sequences[0].seq.size();
  
  StateSeq s(end_sequences[0]);
  if (!FIX_ROOT)
    the_model.sample_state_sequence_init(n_sites, gen, s.seq);
  vector<StateSeq> sequences(n_nodes, s);
  
  bool pass = true;
  for (size_t node_id = 1; pass && (node_id < n_nodes); ++node_id) {
    // start new branch
    vector<GlobalJump> global_path;
    const double curr_branch_len = th.branches[node_id];
    
    TripletSampler ts(sequences[th.parent_ids[node_id]]);
    double time_value = 0;
    while (time_value < curr_branch_len)
      sample_jump(the_model, curr_branch_len, gen, ts, global_path,
                  time_value);
    ts.get_sequence(sequences[node_id]);
    
    // Convert global jumps to local paths
    vector<Path> path_by_site(n_sites);
    for (size_t site_id = 0; site_id < n_sites; site_id++) {
      path_by_site[site_id].init_state = s.seq[site_id];
      path_by_site[site_id].tot_time = curr_branch_len;
    }
    assign_changes_to_sites(global_path, path_by_site);
    paths[node_id] = path_by_site;
    
    // leaf sequence agree?
    pass = th.is_leaf(node_id) ?
    check_leaf_seq(end_sequences[node_id], paths[node_id]) : true;
  }
  return pass;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <class T>
static void
print_vec(vector<T> vec) {
  for (size_t i = 0; i < vec.size(); i++)
    cerr << vec[i] << ", ";
  cerr << endl;
}


/* generate initial paths by heuristics */
static void
initialize_paths(std::mt19937 &gen, const TreeHelper &th,
                 vector<StateSeq> &end_sequences,
                 vector<vector<Path> > &paths, const bool FIX_ROOT) {
  
  const size_t n_sites = end_sequences.front().seq.size();
  
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
          child_states[n_ch++] = (end_sequences[*c].seq[site_id] == '1');
        
        // sample parent state
        if (th.parent_ids[node_id] > 0 || !FIX_ROOT)
          end_sequences[node_id].seq[site_id] = '0' +
          child_states[std::floor(unif()*n_ch)];
        
        // assign site-specific paths above two child nodes
        for (ChildSet c(th.subtree_sizes, node_id); c.good(); ++c) {
          
          const double len = th.branches[*c];
          paths[*c][site_id].tot_time = len;
          paths[*c][site_id].init_state =
            (end_sequences[node_id].seq[site_id] == '1');
          
          const bool parent_state = (end_sequences[node_id].seq[site_id] == '1');
          if (parent_state != paths[*c][site_id].init_state)
            paths[*c][site_id].jumps.push_back(unif() * len);
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*
static void
update_statistics(vector<vector<double> > &J_total,
                  vector<vector<double> > &D_total,
                  const vector<vector<double> > &J,
                  const vector<vector<double> > &D) {
  for (size_t node_id = 0; node_id < J_total.size(); node_id++)
    for (size_t i = 0; i < J_total[node_id].size(); i++) {
      J_total[node_id][i] += J[node_id][i];
      D_total[node_id][i] += D[node_id][i];
    }
}

static void
divide_statistics(vector<vector<double> > &J_total,
                  vector<vector<double> > &D_total, const size_t n) {
  for (size_t node_id = 0; node_id < J_total.size(); node_id++) {
    transform(J_total[node_id].begin(), J_total[node_id].end(),
              J_total[node_id].begin(),
              std::bind(std::divides<double>(), std::placeholders::_1, n));
    transform(D_total[node_id].begin(), D_total[node_id].end(),
              D_total[node_id].begin(),
              std::bind(std::divides<double>(), std::placeholders::_1, n));
  }
}
*/

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

static void
write_statistics_header(const string outfile) {
  std::ofstream out(outfile.c_str());
  out << "ITR\tSAMPLE\tJ_000\tJ_001\tJ_010\tJ_011\tJ_100\tJ_101\tJ_110\tJ_111\t"
  << "D_000\tD_001\tD_010\tD_011\tD_100\tD_101\tD_110\tD_111" << endl;
}


static void
write_statistics(const string outfile, const size_t itr, const size_t sample,
                 const vector<double> &J, const vector<double> &D) {
  std::ofstream out(outfile.c_str(), std::ofstream::app);
  out << itr << "\t" << sample << "\t";
  for (size_t i = 0; i < 8; i++)
    out << J[i] << "\t";
  for (size_t i = 0; i < 7; i++)
    out << D[i] << "\t";
  out << D[7] << endl;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char **argv) {
  try {
    bool VERBOSE = false;
    bool FIX_ROOT = false;
    bool FULL_PATH = false;  /* use SingleSiteSampler (full path) */

    
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
                           "<param>");
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
    opt_parse.add_opt("FULL_PATH", 'F', "use SingleSiteSampler (full path)",
                      false, FULL_PATH);
    opt_parse.add_opt("fixroot", 'R', "fix root states",
                      false, FIX_ROOT);
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

    
    /* (4) GENERATE ONE PATH */
    if (VERBOSE)
      cerr << "[GENERATE A SAMPLE PATH]" << endl;
    vector<vector<Path> > paths; // along multiple branches
    forward_simulation(the_model, th, n_sites, paths, gen);

  
    /* (5) FORWARD SIMULATION */
    if (VERBOSE)
      cerr << "[FORWARD SIMULATION]" << endl;

    // collect sequences
    StateSeq root_seq;
    collect_init_sequences(paths[1], root_seq);
    vector<StateSeq> end_sequences (n_nodes, root_seq);
    for (size_t node_id = 1; node_id < n_nodes; node_id++) {
      collect_end_sequences(paths[node_id], end_sequences[node_id]);
    }
  
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
  
      bool success = end_cond_forward_simulation(the_model, th, end_sequences,
                                                 sampled_paths, gen, FIX_ROOT);
      if (success) {
        vector<vector<double> > J, D;
        get_sufficient_statistics(sampled_paths, J, D);
        for (size_t node_id = 1; node_id < n_nodes; node_id++)
          write_statistics(statfiles_forward[node_id],
                           0, n_forward_samples_collected,
                           J[node_id], D[node_id]);
        n_forward_samples_collected++;
      }
    }

    /* (6) MCMC */
    if (VERBOSE)
      cerr << "[MCMC USING "
      << (FULL_PATH ? "SINGLESITESAMPLER" : "INTERVALSAMPLER") << "]" << endl;

    // mcmc output files
    vector<string> statfiles_mcmc(n_nodes);
    for (size_t node_id = 1; node_id < n_nodes; node_id++) {
      const string fstat = outfile + "." + th.node_names[node_id] + ".mcmc";
      statfiles_mcmc[node_id] = fstat;
      write_statistics_header(fstat);
    }
    
    // distort/randomize paths
    vector<vector<Path> > mcmc_paths;
    initialize_paths(gen, th, end_sequences, mcmc_paths, FIX_ROOT);

    for (size_t site_id = 1; site_id < n_sites - 1; ++site_id)
      assert(mcmc_paths[1][site_id].end_state() == paths[1][site_id].end_state());

    vector<Path> proposed_path;
    for(size_t batch_id = 0; batch_id < n_mcmc_batches; batch_id++) {
      for (size_t sample_id = 0; sample_id < batch; sample_id++) {
        for (size_t site_id = 1; site_id < n_sites - 1; ++site_id) {

          if (FULL_PATH)
            Metropolis_Hastings_site(the_model, th, site_id, mcmc_paths, gen,
                                     proposed_path, FIX_ROOT);
          else
            Metropolis_Hastings_interval(the_model, th, site_id, mcmc_paths,
                                         gen);
          
          assert(mcmc_paths[1][site_id].end_state() == paths[1][site_id].end_state());
        }
        
        // write stats
        vector<vector<double> > J, D;
        get_sufficient_statistics(mcmc_paths, J, D);
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
