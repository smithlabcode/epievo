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
#include <random>
#include <functional>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "Path.hpp"
#include "EpiEvoModel.hpp"
#include "EndCondSampling.hpp"
#include "SingleSiteSampler.hpp"
#include "ContinuousTimeMarkovModel.hpp"
#include "TripletSampler.hpp"
#include "GlobalJump.hpp"
#include "epievo_utils.hpp"


using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::runtime_error;
using std::numeric_limits;
using std::begin;
using std::end;

using std::exponential_distribution;
using std::placeholders::_1;
using std::function;


template <typename T>
using disc_distr = std::discrete_distribution<T>;

void
add_sufficient_statistics(const vector<Path> &paths,
                          vector<double> &J, vector<double> &D) {
  // iterate over sites with valid triples (i.e. not the first and last)
  const size_t n_sites = paths.size();
  for (size_t i = 1; i < n_sites - 1; ++i)
    add_sufficient_statistics(paths[i-1], paths[i], paths[i+1], J, D);
}

void
add_sufficient_statistics(const vector<vector<Path> > &paths,
                          vector<double> &J, vector<double> &D) {
  // iterate over sites with valid triples (i.e. not the first and last)
  const size_t n_sites = paths.size();
  for (size_t i = 1; i < n_sites - 1; ++i)
    add_sufficient_statistics(paths[i-1][1], paths[i][1], paths[i+1][1], J, D);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

static void
collect_init_sequences(const vector<Path> &paths, vector<bool> &seq) {
  seq = vector<bool>(paths.size());
  for (size_t i = 0; i < paths.size(); ++i)
    seq[i] = paths[i].init_state;
}

static void
collect_end_sequences(const vector<Path> &paths, vector<bool> &seq) {
  seq = vector<bool>(paths.size());
  for (size_t i = 0; i < paths.size(); i++)
    seq[i] = paths[i].end_state();
}

/* This function does the sampling for an individual change in the
 * state sequence
 */
static void
sample_jump(const EpiEvoModel &the_model, const double total_time,
            std::mt19937 &gen, TripletSampler &ts,
            vector<GlobalJump> &the_path, double &time_value) {

  static const size_t n_triplets = 8;

  // triplet_count = c_{ijk} for current sequence (encoded in the
  // TripletSampler object)
  vector<size_t> triplet_counts;
  ts.get_triplet_counts(triplet_counts);

  // holding_rate = c_{ijk}*lambda_{ijk}
  const double holding_rate =
    std::inner_product(begin(triplet_counts), end(triplet_counts),
                       begin(the_model.triplet_rates), 0.0);
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
    disc_distr<size_t> multinom(begin(triplet_prob), end(triplet_prob));
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
check_leaf_seq(const vector<bool> &s, const vector<Path> &by_site) {
  bool leaf_is_good = true;
  for (size_t i = 0; i < by_site.size() && leaf_is_good; ++i)
    leaf_is_good = (by_site[i].end_state() == s[i]);
  return leaf_is_good;
}


static void
global_to_local(const vector<bool> &root_seq, const double evo_time,
                const vector<GlobalJump> &gp, vector<Path> &paths) {

  const size_t n_sites = root_seq.size();

  paths = vector<Path>(n_sites);
  for (size_t i = 0; i < n_sites; ++i) {
    paths[i].init_state = root_seq[i];
    paths[i].tot_time = evo_time;
  }

  for (size_t i = 0; i < gp.size(); ++i)
    paths[gp[i].position].jumps.push_back(gp[i].timepoint);
}


static void
forward_simulation(const EpiEvoModel &the_model,
                   const double &curr_branch_len,
                   const size_t n_sites, vector<Path> &paths,
                   std::mt19937 &gen) {

  // sample a root
  vector<bool> root_seq;
  the_model.sample_state_sequence_init(n_sites, gen, root_seq);
  TripletSampler ts(root_seq);

  // generate the global path
  vector<GlobalJump> global_path;
  double time_value = 0;
  while (time_value < curr_branch_len)
    sample_jump(the_model, curr_branch_len, gen, ts, global_path, time_value);

  // convert global to site-specific histories
  global_to_local(root_seq, curr_branch_len, global_path, paths);
}


static bool
end_cond_forward_simulation(const EpiEvoModel &the_model, const double &evo_time,
                            const vector<bool> &root_seq,
                            const vector<bool> &leaf_seq,
                            vector<Path> &paths, std::mt19937 &gen) {

  // start new branch
  const double curr_branch_len = evo_time;
  TripletSampler ts(root_seq);

  vector<GlobalJump> global_path;
  double time_value = 0.0;
  while (time_value < curr_branch_len)
    sample_jump(the_model, curr_branch_len, gen, ts, global_path, time_value);

  vector<bool> s;
  ts.get_sequence(s);

  // Convert global jumps to local paths
  global_to_local(root_seq, curr_branch_len, global_path, paths);

  // leaf sequence agree?
  const bool leaf_is_good = check_leaf_seq(leaf_seq, paths);

  return leaf_is_good;
}


/* generate initial paths by heuristics */
static void
initialize_paths(const vector<bool> &root_seq, const vector<bool> &leaf_seq,
                 const double evo_time, vector<vector<Path> > paths,
                 std::mt19937 &gen) {
  
  paths.resize(root_seq.size());
  
  auto unif =
  bind(std::uniform_real_distribution<double>(0.0, 1.0), std::ref(gen));
  
  for (size_t site_id = 0; site_id < root_seq.size(); ++site_id) {
    paths[site_id].resize(2);
    paths[site_id][1].init_state = root_seq[site_id];
    paths[site_id][1].tot_time = evo_time;
    
    if (root_seq[site_id] != leaf_seq[site_id])
      paths[site_id][1].jumps.push_back(unif() * evo_time);
  }
}

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
  std::ofstream out(outfile, std::ofstream::app);
  out << itr << "\t" << sample << "\t";
  copy(begin(J), begin(J) + 8, std::ostream_iterator<double>(out, "\t"));
  copy(begin(D), begin(D) + 8, std::ostream_iterator<double>(out, "\t"));
  out << endl;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char **argv) {
  try {
    bool VERBOSE = false;

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
    opt_parse.add_opt("evo-time", 'T', "evolutionary time", true,
                      evolutionary_time);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);;
    opt_parse.add_opt("outfile", 'o', "outfile (prefix)",
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
    ///////////////////////////////////////////////////////////////////////////

    /* (1) READING PARAMETERS FROM FILE */
    if (VERBOSE)
      cerr << "[READING PARAMETERS: " << param_file << "]" << endl;
    EpiEvoModel the_model;
    read_model(param_file, the_model);
    the_model.scale_triplet_rates();
    if (VERBOSE)
      cerr << the_model << endl;

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
      cerr << "[GENERATE A SAMPLE PATH]"<< endl;
    vector<Path> paths; // paths on single branch
    forward_simulation(the_model, evolutionary_time, n_sites, paths, gen);
    // collect sequences
    vector<bool> root_seq;
    collect_init_sequences(paths, root_seq);
    
    vector<bool> leaf_seq;
    collect_end_sequences(paths, leaf_seq);
    if (VERBOSE) {
      cout << "ROOT SEQ: ";
      copy(begin(root_seq), end(root_seq),
           std::ostream_iterator<int>(cout, ""));
      cout << "\nEND SEQ: ";
      copy(begin(leaf_seq), end(leaf_seq),
           std::ostream_iterator<int>(cout, ""));
      cout << endl;
    }
    
    /* (5) FORWARD SIMULATION */
    if (VERBOSE)
      cerr << "[FORWARD SIMULATION]" << endl;
    // forward-simulation output files
    const string fstat = outfile + ".forward";
    write_statistics_header(fstat);

    // sampling
    cerr << "[SAMPLING]" << endl;
    size_t n_forward_samples_collected = 0;
    while (n_forward_samples_collected < batch) {

      vector<Path> sampled_paths;
      bool success = end_cond_forward_simulation(the_model, evolutionary_time,
                                                 root_seq, leaf_seq,
                                                 sampled_paths, gen);
      if (success) {
        vector<double> J(8), D(8);
        add_sufficient_statistics(sampled_paths, J, D);
        write_statistics(fstat, 0, n_forward_samples_collected, J, D);
        n_forward_samples_collected++;
      }
    }

    if (VERBOSE)
      cerr << "[MCMC USING SINGLESITESAMPLER]" << endl;

    // mcmc output files
    const string mcmc_stat = outfile + ".mcmc";
    write_statistics_header(mcmc_stat);

    // distort/randomize paths
    vector<vector<Path> > mcmc_paths(2); // [sites] x [nodes]
    initialize_paths(root_seq, leaf_seq, evolutionary_time, mcmc_paths, gen);
    const TreeHelper th(evolutionary_time);

    // [sites] x [nodes] x [emit states]
    vector<vector<vector<double> > > emit(n_sites);
    for (size_t site_id = 0; site_id < n_sites; site_id++) {
      emit[site_id].resize(2);
      emit[site_id][0].resize(2);
      emit[site_id][0][0] = (mcmc_paths[1][site_id].init_state == false ?
                             1.0 : 0.0);
      emit[site_id][0][1] = (mcmc_paths[1][site_id].init_state == true ?
                             1.0 : 0.0);
    }

    vector<double> tri_llh(n_sites, 0.0);
    // pre-compute triplet log-likelihood on each site
    for (size_t site_id = 1; site_id < n_sites - 1; ++site_id)
      tri_llh[site_id] = path_log_likelihood(the_model, mcmc_paths[site_id-1],
                                             mcmc_paths[site_id],
                                             mcmc_paths[site_id+1]);

    for (size_t batch_id = 0; batch_id < n_mcmc_batches; batch_id++) {
      for (size_t sample_id = 0; sample_id < batch; sample_id++) {
        for (size_t site_id = 1; site_id < n_sites - 1; ++site_id) {
          Metropolis_Hastings_site(the_model, th, site_id, mcmc_paths,
                                   emit[site_id], tri_llh[site_id-1],
                                   tri_llh[site_id], tri_llh[site_id+1], gen);
        }

        // write stats
        vector<double> J(8), D(8);
        add_sufficient_statistics(mcmc_paths, J, D);
        write_statistics(mcmc_stat, batch_id, sample_id, J, D);
      }
    }

  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
