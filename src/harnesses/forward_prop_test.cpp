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
#include <random>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <istream>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "PhyloTreePreorder.hpp"
#include "TreeHelper.hpp"
#include "TripletSampler.hpp"
#include "StateSeq.hpp"
#include "EpiEvoModel.hpp"
#include "GlobalJump.hpp"
#include "Path.hpp"
#include "ParamEstimation.hpp"

using std::vector;
using std::string;
using std::endl;
using std::cerr;
using std::cout;
using std::to_string;
using std::ostream_iterator;
using std::numeric_limits;
using std::istringstream;
using std::runtime_error;


////////////////////////////////////////////////////////////////////////////////
static void
delete_last_branch(vector<vector<Path> > &paths, vector<string> &node_names,
                   TreeHelper &th) {
  paths.pop_back();
  node_names.pop_back();
  th.branches.pop_back();
  th.parent_ids.pop_back();
  th.subtree_sizes.pop_back();
  th.subtree_sizes[0]--;
  th.n_nodes--;
}

bool
file_is_readable(const string &param_file) {
  std::ifstream in(param_file.c_str());
  return in.good();
}


static void
write_root_to_pathfile_local(const string &outfile, const string &root_name) {
  std::ofstream outpath(outfile.c_str());
  outpath << "NODE:" << root_name << endl;
}


static void
append_to_pathfile_local(const string &pathfile, const string &node_name,
                         const vector<Path> &path_by_site) {
  std::ofstream outpath(pathfile.c_str(), std::ofstream::app);
  outpath << "NODE:" << node_name << endl;
  for (size_t i = 0; i < path_by_site.size(); ++i)
    outpath << i << '\t' << path_by_site[i] << '\n';
}


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

static bool
check_leaf_seq(const StateSeq s, const vector<Path> &by_site) {
  bool pass = true;
  size_t site_id = 0;
  while(pass && site_id < by_site.size()) {
    const bool leaf_state = s.seq[site_id];
  //   cerr << leaf_state << endl;
    pass = (by_site[site_id].end_state() == leaf_state) ? true : false;
    site_id++;
  }
  return pass;
}

static bool
check_homo_sites(const vector<Path> &by_site,
                 const vector<size_t> censored_site) {
  bool homo = true;
  size_t site_id = 0;
  while(homo && site_id < by_site.size()) {
    if (std::find(censored_site.begin(), censored_site.end(),
                  site_id) == censored_site.end()) {
      homo = (by_site[site_id].jumps.empty()) ? true : false;
    }
    site_id++;
  }
  return homo;
}


int main(int argc, const char **argv) {

  try {
    
    bool VERBOSE = false;
    bool FIX_NEIGHBORS = false;
    string outfile;
    string statfile;
    
    size_t rounds = 10;
    size_t rng_seed = std::numeric_limits<size_t>::max();
    size_t the_branch = 1;
    vector<size_t> consored_site = {2};

    ////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]),
                           "Forward posterior sampling"
                           "from (stationary) simulation",
                           "<params-file>");
    opt_parse.add_opt("rounds", 'r', "number of MCMC rounds",
                      false, rounds);
    opt_parse.add_opt("branch", 'b', "branch to simulate",
                      false, the_branch);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);
    opt_parse.add_opt("statfile", 'S', "MCMC stats",
                      false, statfile);
    opt_parse.add_opt("outfile", 'o', "outfile (sampling paths)",
                      false, outfile);
    opt_parse.add_opt("fix", 'f', "fix neighbors",
                      false, FIX_NEIGHBORS);
    
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
    
    const string param_file(leftover_args[0]);
    const string treefile(leftover_args[1]);
    const string input_file(leftover_args[2]);
    ////////////////////////////////////////////////////////////////////////
    /* (1) LOADING (FAKE) TREE */
    if (VERBOSE)
      cerr << "[READING TREE: " << treefile << "]" << endl;
    PhyloTreePreorder the_tree; // tree topology and branch lengths
    std::ifstream tree_in(treefile.c_str());
    if (!tree_in || !(tree_in >> the_tree))
      throw runtime_error("cannot read tree file: " + treefile);
    TreeHelper th(the_tree);

    /* (2) READING PARAMETERS FROM FILE */
    if (VERBOSE)
      cerr << "[READING PARAMETERS: " << param_file << endl;
    EpiEvoModel the_model;
    read_model(param_file, the_model);
    the_model.scale_triplet_rates();
    if (VERBOSE)
      cerr << the_model << endl;

    /* (3) LOADING PATHS */
    if (VERBOSE)
      cerr << "[READING PATHS FILE: " << input_file << "]" << endl;
    vector<string> node_names;
    vector<vector<Path> > paths; // along multiple branches
    read_paths(input_file, node_names, paths);
    const size_t n_sites = paths[the_branch].size();
    // remove unwanted branch for now
    delete_last_branch(paths, node_names, th);
    size_t n_nodes = paths.size();
    
    vector<char> root_seq;
    if (VERBOSE)
      cerr << "[INPUT ROOT SEQUENCE: ]" << endl;
    for (size_t pos = 0; pos < paths[the_branch].size(); pos++) {
      cerr << paths[the_branch][pos].init_state;
      root_seq.push_back(paths[the_branch][pos].init_state);
    }
    cerr << endl;
  
    if (VERBOSE)
      cerr << "[INPUT LEAF SEQUENCE: ]" << endl;
    for (size_t pos = 0; pos < paths[the_branch].size(); pos++) {
      cerr << paths[the_branch][pos].end_state();
    }
    cerr << endl;
  
    /* (4) INITIALIZING THE RANDOM NUMBER GENERATOR */
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    std::mt19937 gen(rng_seed);
    
    /* (5) POSTERIOR SAMPLING */
    string rootfile = statfile + ".root";
    std::ofstream of_root;
    of_root.open(rootfile.c_str());
    std::ofstream outstat(statfile.c_str());
    std::ostream out_root(of_root.rdbuf());
    outstat << "J_000\tJ_001\tJ_010\tJ_011\tJ_100\tJ_101\tJ_110\tJ_111\t"
    << "D_000\tD_001\tD_010\tD_011\tD_100\tD_101\tD_110\tD_111" << endl;
    
    size_t num_sample_collected = 0;
    size_t num_single_process = 0;
    while (num_sample_collected < rounds) {
      bool accept = false;
      vector<vector<Path> > updated_paths(paths);

      /* (a) SAMPLING ROOT SEQUENCE */
      // vector<char> root_seq;
      the_model.sample_state_sequence_init(n_sites, gen, root_seq);
      
      StateSeq s(root_seq);
      vector<StateSeq> sequences(n_nodes, s);
      vector<vector<GlobalJump> > global_paths;

      for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
        vector<GlobalJump> global_path;
        const double curr_branch_len = th.branches[node_id];
   

        TripletSampler ts(sequences[th.parent_ids[node_id]]);
        double time_value = 0;
        while (time_value < curr_branch_len)
          sample_jump(the_model, curr_branch_len, gen, ts, global_path,
                      time_value);

        ts.get_sequence(sequences[node_id]);
        num_single_process++;
        
        // Convert global jumps to local paths
        vector<Path> path_by_site(n_sites);
        for (size_t site_id = 0; site_id < n_sites; site_id++) {
          path_by_site[site_id].init_state = s.seq[site_id];
          path_by_site[site_id].tot_time = curr_branch_len;
        }
        assign_changes_to_sites(global_path, path_by_site);
        updated_paths[node_id] = path_by_site;
        
        // leaf sequence agree?
        if (check_leaf_seq(sequences[node_id], paths[node_id]) &&
            (!FIX_NEIGHBORS || check_homo_sites(path_by_site, consored_site))) {
          accept = true;
          global_paths.push_back(global_path);
          
          // output
          vector<double> J;
          vector<double> D;
          get_sufficient_statistics(updated_paths, J, D);
          for (size_t i = 0; i < 8; i++)
            outstat << J[i] << "\t";
          for (size_t i = 0; i < 7; i++)
            outstat << D[i] << "\t";
          outstat << D[7] << endl;
          // write root
          for (size_t site_id = 0; site_id < n_sites - 1; site_id++)
            out_root << updated_paths[the_branch][site_id].init_state << "\t";
          out_root << updated_paths[the_branch][n_sites - 1].init_state << endl;

          num_sample_collected++;
        }
      }
    }
    
    cerr << "SUCCESS RATE: "
    << static_cast<double> (num_sample_collected) / num_single_process << endl;
    
    /* (6) OUTPUT */
    write_root_to_pathfile_local(outfile, th.node_names.front());
    for (size_t node_id = 1; node_id < th.n_nodes; ++node_id)
      append_to_pathfile_local(outfile, th.node_names[node_id], paths[node_id]);
  }
  
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
