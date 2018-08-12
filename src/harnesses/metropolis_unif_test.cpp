/* Copyright (C) 2018 University of Southern California
 *                    Liz Ji, Jianghan Qu and Andrew D Smith
 *
 * single_site_sampling_test: This program is to test the
 * functionality in SingleSiteSampler.*pp which includes a combination
 * of end-point sampling, in an up-down process, with virtual nodes
 * internal to a branch, along with end-point conditioned sampling.
 *
 * Author: Liz Ji, Andrew D. Smith and Jianghan Qu
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
#include "TripletSampler.hpp"
#include "GlobalJump.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::runtime_error;
using std::numeric_limits;


static void
assign_changes_to_sites(const vector<GlobalJump> &global_path,
                        vector<Path> &by_site) {

  const size_t n_changes = global_path.size();
  for (size_t i = 0; i < n_changes; ++i)
    by_site[global_path[i].position].jumps.push_back(global_path[i].timepoint);
}

///////////////////////////////////////////////////////////////////////////

/* This function does the sampling for an individual change in the
   state sequence
 */
static size_t
count_jumps(const vector<vector<Path> > paths, const size_t the_site) {
  size_t counts = 0;
  for (size_t b = 1; b < paths.size(); ++b)
    counts += paths[b][the_site].jumps.size();
  return counts;
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


int main(int argc, const char **argv) {
  try {
    bool VERBOSE = false;
    string outfile;
    string statfile;

    size_t rounds = 10;
    size_t rng_seed = std::numeric_limits<size_t>::max();
    size_t the_branch = 1;
    ///////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]), "test mcmc procedure",
                           "<param> <treefile> <path_file>");
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
    const string treefile(leftover_args[1]);
    const string input_file(leftover_args[2]);
    ///////////////////////////////////////////////////////////////////////////
    /* (1) LOADING (FAKE) TREE */
    if (VERBOSE)
      cerr << "[READING TREE: " << treefile << "]" << endl;
    PhyloTreePreorder the_tree; // tree topology and branch lengths
    std::ifstream tree_in(treefile.c_str());
    if (!tree_in || !(tree_in >> the_tree))
      throw runtime_error("cannot read tree file: " + treefile);
    const TreeHelper th(the_tree);

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
    
    /* (4) INITIALIZING THE RANDOM NUMBER GENERATOR */
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    std::mt19937 gen(rng_seed);
    if (VERBOSE)
      cerr << "[INITIALIZED RNG (seed=" << rng_seed << ")]" << endl;

    /* (5) MCMC */
    std::ofstream outstat(statfile.c_str());
    outstat << "ROUND\tSITE\tNJUMPS_OLD\tNJUMPS_NEW\tACCEPTED" << endl;
    for (size_t r = 1; r <= rounds; ++r) {
      if (VERBOSE)
        cerr << "*********************\nROUND: " << r << endl;
      for (size_t site_id = 2; site_id < n_sites - 2; ++site_id) {
        outstat << r << "\t" << site_id << "\t"
                << count_jumps(paths, site_id) << "\t";
        if (VERBOSE)
          cerr << "site: " << site_id
               << ", original n_jumps:" << count_jumps(paths, site_id);
        vector<Path> proposed_path;
        const bool accpeted = Metropolis_Hastings_site(the_model, th, site_id,
                                                       paths, gen,
                                                       proposed_path, r == 1);
        outstat << count_jumps(paths, site_id) << "\t" << accpeted << endl;
        
        if (VERBOSE)
          cerr << ", updated n_jumps:" << count_jumps(paths, site_id) << endl;
        
      }
    }
    
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
