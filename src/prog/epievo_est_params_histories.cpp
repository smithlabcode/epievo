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
#include <cmath> /* log */

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "Path.hpp"
#include "EpiEvoModel.hpp"
#include "EndCondSampling.hpp"
#include "ContinuousTimeMarkovModel.hpp"
#include "TreeHelper.hpp"
#include "SingleSiteSampler.hpp"
#include "TripletSampler.hpp"
#include "GlobalJump.hpp"
#include "ParamEstimation.hpp"
#include "epievo_utils.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::runtime_error;
using std::numeric_limits;
using std::ofstream;

static void
write_root_to_pathfile_local(const string &outfile, const string &root_name) {
  ofstream out(outfile);
  if (!out)
    throw std::runtime_error("bad output file: " + outfile);
  out << "NODE:" << root_name << endl;
}

static void
append_to_pathfile_local(const string &pathfile, const string &node_name,
                         const vector<vector<Path> > &paths,
                         const size_t node_id) {
  std::ofstream out(pathfile, std::ofstream::app);
  if (!out)
    throw std::runtime_error("bad output file: " + pathfile);

  out << "NODE:" << node_name << endl;
  for (size_t i = 0; i < paths.size(); ++i)
    out << i << '\t' << paths[i][node_id] << '\n';
}

static void
write_mcmc_verbose_header(std::ostream &out) {
  vector<string> header_tokens = {
    "itr",
    "stationary",
    "baseline",
    "init",
    "acc_rate",
    "llh",
  };
  copy(begin(header_tokens), end(header_tokens),
       std::ostream_iterator<string>(out, "\t"));
  out << '\n';
}


int main(int argc, const char **argv) {
  try {

    bool VERBOSE = false;
    bool optimize_branches = false;
    bool ONEBRANCH = false;

    string outfile;
    string param_file_updated;
    string treefile_updated;

    size_t iteration = 10;           // MCMC-EM iterations
    size_t batch = 10;              // MCMC iterations
    size_t burnin = 10;           // burn-in MCMC iterations
    double increment_k = 0.0;           // batch increment controller

    size_t rng_seed = std::numeric_limits<size_t>::max();

    static const double param_tol = 1e-10;
    ///////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]),
                           "estimate parameters and evolutionary histories",
                           "<param> <treefile> <path_file>");
    opt_parse.add_opt("iteration", 'i', "number of MCMC-EM iteration",
                      false, iteration);
    opt_parse.add_opt("batch", 'B', "number of MCMC iteration", false, batch);
    opt_parse.add_opt("increment", 'k', "batch increment controller",
                      false, increment_k);
    opt_parse.add_opt("burnin", 'L', "MCMC burn-in length",
                      false, burnin);
    opt_parse.add_opt("one-branch", 'T', "one-branch tree", false,
                      ONEBRANCH);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);
    opt_parse.add_opt("outfile", 'o', "output file of local paths",
                      true, outfile);
    opt_parse.add_opt("outparam", 'p', "output file of parameters",
                      false, param_file_updated);
    opt_parse.add_opt("outtree", 't', "output file of tree",
                      false, treefile_updated);
    opt_parse.add_opt("branch", 'b', "optimize branch lengths",
                      false, optimize_branches);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.set_show_defaults();
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
    if (leftover_args.size() < 3) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    const string param_file(leftover_args[0]);
    const string tree(leftover_args[1]);
    const string input_file(leftover_args[2]);

    ///////////////////////////////////////////////////////////////////////////
    /* (1) LOADING (FAKE) TREE */
    if (VERBOSE)
      cerr << "[READING TREE: " << tree << "]" << endl;
    PhyloTreePreorder the_tree; // tree topology and branch lengths
    TreeHelper th;
    if (ONEBRANCH) {
      if (VERBOSE)
        cerr << "initializing two node tree with time: " << tree << endl;
      th = TreeHelper(std::stod(tree));
    }
    else {
      cerr << "reading tree file: " << tree << endl;
      std::ifstream tree_in(tree);
      if (!tree_in || !(tree_in >> the_tree))
        throw std::runtime_error("bad tree file: " + tree);
      th = TreeHelper(the_tree);
    }

    /* (2) READING PARAMETERS FROM FILE */
    if (VERBOSE)
      cerr << "[READING PARAMETERS: " << param_file << "]" << endl;
    EpiEvoModel the_model;
    read_model(param_file, the_model);
    if (VERBOSE)
      cerr << the_model << endl;

    /* (3) LOADING PATHS */
    if (VERBOSE)
      cerr << "[READING PATHS FILE: " << input_file << "]" << endl;
    vector<string> node_names;
    vector<vector<Path> > orig_paths; // [nodes] x [sites]
    vector<vector<Path> > paths; // [sites] x [nodes]
    read_paths(input_file, node_names, orig_paths);

    ////////// transposing
    const size_t n_sites = orig_paths[1].size();
    const size_t n_nodes = orig_paths.size();
    paths.resize(n_sites);
    for (size_t i = 0; i < n_sites; ++i) {
      paths[i].resize(n_nodes);
      for (size_t b = 1; b < n_nodes; ++b) {
        paths[i][b] = orig_paths[b][i];
      }
    }
    //////////

    /* (4) INITIALIZING THE RANDOM NUMBER GENERATOR */
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    if (VERBOSE)
      cerr << "rng seed: " << rng_seed << endl;
    std::mt19937 gen(rng_seed);

    /* (5) MCMC */

    /* MCMC PARAMETERS AND DATA*/
    vector<vector<double> > J;
    vector<vector<double> > D;
    two_by_two root_frequencies;
    vector<double> tri_llh(n_sites, 0.0);

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    /////// START MCMC-EM
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    SingleSiteSampler mcmc(th.n_nodes);

    if (VERBOSE) {
      write_mcmc_verbose_header(cerr);
      cerr << "0" << "\t" << the_model.T(0, 0) << "\t"
      << the_model.T(1, 1) << "\t"
      << the_model.stationary_baseline(0, 0) << "\t"
      << the_model.stationary_baseline(1, 1) << endl;
    }

    size_t current_batch = batch;

    /* METROPOLIS-HASTINGS ALGORITHM */
    for (size_t itr = 0; itr  < iteration; itr++) {

      // pre-compute log(rates)
      double log_rates[8];
      std::transform(std::begin(the_model.triplet_rates),
                     std::end(the_model.triplet_rates),
                     std::begin(log_rates),
                     static_cast<double(*)(double)>(log));

      // pre-compute triplet log-likelihood on each site
      for (size_t site_id = 1; site_id < n_sites - 1; ++site_id)
        tri_llh[site_id] = path_log_likelihood(the_model, paths[site_id-1],
                                               paths[site_id], paths[site_id+1],
                                               log_rates);

      // Burning
      for(size_t burnin_itr = 0; burnin_itr < burnin; burnin_itr++) {
        for (size_t site_id = 1; site_id < n_sites - 1; ++site_id) {
          mcmc.Metropolis_Hastings_site(the_model, th, site_id, paths,
                                        tri_llh[site_id-1], tri_llh[site_id],
                                        tri_llh[site_id+1], gen, log_rates);
        }
      }

      // MCMC samples
      current_batch += (itr + 1)*increment_k;
      size_t n_accepted = 0;
      vector<vector<double> > J_accum(th.n_nodes, vector<double>(8, 0.0));
      vector<vector<double> > D_accum(th.n_nodes, vector<double>(8, 0.0));

      two_by_two root_accum;

      for (size_t mcmc_itr = 0; mcmc_itr < current_batch; mcmc_itr++) {
        for (size_t site_id = 1; site_id < n_sites - 1; ++site_id)
          n_accepted += mcmc.Metropolis_Hastings_site(the_model, th, site_id,
                                                      paths,
                                                      tri_llh[site_id-1],
                                                      tri_llh[site_id],
                                                      tri_llh[site_id+1],
                                                      gen, log_rates);
        /* CALCULATE SUFFICIENT STATS */
        get_sufficient_statistics(paths, J, D);
        get_root_frequencies(paths, root_frequencies);
        /* RECORD STATS (POST-BURN-IN) */
        for (size_t b = 1; b < J.size(); b++) {
          for (size_t i = 0; i < 8; i++) {
            J_accum[b][i] += J[b][i];
            D_accum[b][i] += D[b][i];
          }
        }
        root_accum += root_frequencies;
      }

      //////////////////////////////////////////////////////////////////////////
      ///////////               POST-MCMC PROCESSING                 ///////////
      //////////////////////////////////////////////////////////////////////////
      /* CALCULATE/OUTPUT BATCH AVERAGE */
      for (size_t b = 1; b < J.size(); ++b) {
        for (size_t i = 0; i < 8; i++) {
          J_accum[b][i] /= current_batch;
          D_accum[b][i] /= current_batch;
        }
      }
      root_accum.divide_all(current_batch);

      /* PARAMETER ESTIMATION */
      double llh = 0.0;
      if (!optimize_branches)
        llh = estimate_rates(false, param_tol, J_accum, D_accum, the_model);
      else {
        llh = estimate_rates_and_branches(false, param_tol, J_accum,
                                          D_accum, th, the_model);
        scale_jump_times(paths, th);
        the_tree.set_branch_lengths(th.branches);
      }

      if (VERBOSE)
        cerr << itr+1 << "\t" << the_model.T(0, 0) << "\t"
             << the_model.T(1, 1) << "\t"
             << the_model.stationary_baseline(0, 0) << "\t"
             << the_model.stationary_baseline(1, 1) << "\t"
             << static_cast<double>(n_accepted)/(current_batch*(n_sites - 2))
             << "\t" << llh << endl;

      /* WRITE PARAMETER FILE */
      if (!param_file_updated.empty()) {
        ofstream out_param(param_file_updated);
        if (!out_param)
          throw runtime_error("bad output param file: " + param_file_updated);
        out_param << the_model.format_for_param_file() << endl;
      }

      /* WRITE PATH FILE */
      write_root_to_pathfile_local(outfile, th.node_names.front());
      for (size_t node_id = 1; node_id < th.n_nodes; ++node_id)
        append_to_pathfile_local(outfile, th.node_names[node_id], paths, node_id);

      /* WRITE TREE FILE */
      if (optimize_branches)
        if (!treefile_updated.empty()) {
          ofstream out_tree(treefile_updated);
          if (!out_tree)
            throw runtime_error("bad output param file: " + treefile_updated);
          out_tree << the_tree << endl;
        }
    }
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
