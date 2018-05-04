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

////////////////////////////////////////////////////////////////////////////////
//////////   copied from global_jumps_to_paths.cpp                    //////////
////////////////////////////////////////////////////////////////////////////////
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


////////////////////////////////////////////////////////////////////////////////
void
paths_to_states_start(const vector<Path> &path_by_site,
                      vector<size_t> &states) {
  states.clear();
  for (size_t i = 0; i < path_by_site.size(); ++i) {
    states.push_back(path_by_site[i].init_state ? 1 : 0);
  }
  assert(states.size() == path_by_site.size());
}

void
paths_to_states_end(const vector<Path> &path_by_site,
                    vector<size_t> &states) {
  states.clear();
  for (size_t i = 0; i < path_by_site.size(); ++i) {
    states.push_back((size_t)(path_by_site[i].end_state()));
  }
}

static void
write_states_at_nodes(const string &outfile, const EpiEvoModel &the_model,
                      const vector<string> &node_names,
                      const vector<vector<Path> > &all_paths) {
  std::ofstream out(outfile.c_str());
  if (!out)
    throw std::runtime_error("bad output file: " + outfile);

  const size_t n_sequences = all_paths.size();
  const size_t n_sites = all_paths[1].size();
  cerr << n_sequences << " nodes and " << n_sites << " sites" << endl;
  out << '#';
  copy(node_names.begin(), node_names.end(),
       std::ostream_iterator<string>(out, "\t"));
  out << '\n';

  vector<vector<size_t> > states(n_sequences);
  paths_to_states_start(all_paths[1], states[0]);
  for (size_t i = 1; i < n_sequences; ++i) {
    paths_to_states_end(all_paths[i], states[i]);
  }

  for (size_t i = 0; i < n_sites; ++i) {
    out << i;
    for (size_t j = 0; j < n_sequences; ++j)
      out << '\t' << states[j][i] ;
    out << '\n';
  }

}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char **argv) {
  try {

    bool VERBOSE = false;
    bool SCALE = true;
    string outfile;
    string outstatefile;

    size_t rounds = 1;

    string param_file;
    string tree_file;

    size_t rng_seed = std::numeric_limits<size_t>::max();

    ///////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]), "test triple path",
                           " <paths-file>");
    opt_parse.add_opt("param", 'p', "parameter file",
                      true, param_file);
    opt_parse.add_opt("tree", 't', "tree file in newick format",
                      true, tree_file);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);
    opt_parse.add_opt("rounds", 'r', "number of posterior update cycles",
                      true, rounds);
    opt_parse.add_opt("outfile", 'o', "outfile for posterior-updated paths",
                      true, outfile);
    opt_parse.add_opt("outstatefile", 'S', "outfile for states at nodes",
                      false, outstatefile);
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string pathsfile(leftover_args.front());
    ///////////////////////////////////////////////////////////////////////////

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
      cerr << "[READING PARAMETER FILE: " << param_file << ", "
           << tree_file << "]" << endl;

    EpiEvoModel the_model;
    read_model(SCALE, param_file, tree_file, the_model);

    if (VERBOSE)
      cerr << the_model << endl;

    /* standard mersenne_twister_engine seeded with rd() */
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }

    if (VERBOSE)
      cerr << "rng seed: " << rng_seed << endl;

    std::mt19937 gen(rng_seed);

    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    cerr << "-- TEST homogeneous interval end-cond. sampling BELOW --" << endl;
    double T = 0.5;
    size_t a = 1;
    size_t b = 1;
    cerr << "Setup: interval of time\t " << T << endl
         << "       neighbor context\t 0-1" << endl
         << "       start state     \t " << a << endl
         << "       end state       \t " << b << endl;

    // prepare helper values
    vector<double> rates(2, 0.0);
    rates[0] = the_model.triplet_rates[0]; // 000
    rates[1] = the_model.triplet_rates[2]; // 010
    vector<double> eigen_vals;
    vector<vector<double> > U;
    vector<vector<double> > Uinv;

    decompose(rates, eigen_vals, U, Uinv);

    vector<double> jump_times;
    end_cond_sample(rates, eigen_vals, U, Uinv,
                    a, b, T, gen, jump_times);

    cerr << " sampled jump times: " << endl;
    for (size_t i = 0; i < jump_times.size(); ++i)
      cerr << jump_times[i] << ",\t" ;
    cerr << endl << "total " << jump_times.size()-2 << " jumps" << endl;

    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    cerr << "----- TEST upward downward sampling BELOW ---------" << endl;

    /*
      size_t site = 10;
      write_root_to_pathfile_local(outfile, the_model.node_names.front());
      for (size_t i = 0; i < 10000000; ++i) {
      vector<Path> new_path;
      gibbs_site(the_model, site, all_paths, gen, new_path);

      std::ofstream outpath(outfile.c_str(), std::ofstream::app);
      for (size_t node_id = 4; node_id < n_nodes; ++node_id) {
      if (new_path[node_id].jumps.size() > 0)
      outpath << "NODE:" << the_model.node_names[node_id] << ";" << i << '\t'
      << site << '\t' << new_path[node_id] << '\n';
      }
      }
    */


    for (size_t k = 0; k < rounds; ++k) {
      cerr << "pass " << k+1 << endl;
      for (size_t i = 2; i < n_sites - 2; ++i) {
        std::uniform_int_distribution<> dis(2, n_sites - 3);
        size_t site = dis(gen);
        vector<Path> new_path;
        gibbs_site(the_model, site, all_paths, gen, new_path);
        for (size_t node_id = 1; node_id < n_nodes; ++node_id)
          all_paths[node_id][site] = new_path[node_id];
      }

      if (!outstatefile.empty()) {
        const string m_outstatefile = outstatefile + std::to_string(k);
        write_states_at_nodes(m_outstatefile, the_model,
                              node_names, all_paths);
      }
    }
    cerr << endl;

    write_root_to_pathfile_local(outfile, the_model.node_names.front());

    for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
      append_to_pathfile_local(outfile, the_model.node_names[node_id],
                               all_paths[node_id]);
    }

  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
