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
#include <iomanip>
#include <sstream>
#include <bitset>

#include <random>
#include <cmath>   /* exp, sqrt, pow */
#include <numeric>  /* std::inner_product */

#include <sys/stat.h>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "PhyloTreePreorder.hpp"
#include "TripletSampler.hpp"
#include "StateSeq.hpp"
#include "EpiEvoModel.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::to_string;


bool
file_exist(const string &param_file) {
  struct stat buf;
  return (stat(param_file.c_str(), &buf) == 0);
}

bool
file_is_readable(const string &param_file) {
  std::ifstream in(param_file.c_str());
  return in.good();
}


static void
count_triplets(const vector<char> &state_sequence,
               vector<size_t> &triplet_count) {
  static const size_t n_triplets = 8;

  triplet_count = vector<size_t>(n_triplets, 0);
  for (size_t i = 1; i < state_sequence.size()-1; ++i)
    ++triplet_count[triple2idx(state_sequence[i-1],
                               state_sequence[i],
                               state_sequence[i+1])];
}


struct Segment {
  Segment(const double tp, const size_t pos) : time_point(tp), position(pos) {}
  double time_point;
  size_t position;
  bool operator<<(const Segment &other) const {
    return time_point < other.time_point;
  }
};

std::ostream &
operator<<(std::ostream &os, const Segment &s) {
  return os << s.time_point << '\t' << s.position;
}

// use new data structure PatSeq
static void
sample_jump(const EpiEvoModel &the_model, const double total_time,
            std::mt19937 &gen, TripletSampler &ts, vector<Segment> &the_path,
            double &time_value) {

  static const size_t n_triplets = 8;

  // triplet_count = c_{ijk} for current sequence (encoded in the
  // TripletSapler object)
  vector<size_t> triplet_counts;
  ts.get_triplet_counts(triplet_counts);

  // holding_rate = c_{ijk}*lambda_{ijk}
  const double holding_rate =
    std::inner_product(triplet_counts.begin(), triplet_counts.end(),
                       the_model.triplet_rates.begin(), 0.0);

  // sample a holding time = time until next state change
  std::exponential_distribution<double> exp_distr(holding_rate);
  const double holding_time = exp_distr(gen);

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
    the_path.push_back(Segment(time_value, change_position));
  }

}


static
void write_pathfile_header(const string &pathfile) {
  std::ofstream outpath(pathfile.c_str());
  outpath << "## paths" << endl;
}

static void
append_to_pathfile(const string &pathfile, const size_t node_id,
                   const vector<Segment> &the_path) {
  std::ofstream outpath(pathfile.c_str(), std::ofstream::app);
  outpath << "NODE\t" << node_id << '\n';
  for (size_t i = 0; i < the_path.size(); ++i)
    outpath << the_path[i] << '\n';
}


static void
write_output(const string &outfile,
             const vector<string> &node_names,
             const vector<vector<char> > &sequences) {

  std::ofstream out(outfile.c_str());
  if (!out)
    throw std::runtime_error("bad output file: " + outfile);

  const size_t n_sequences = sequences.size();
  const size_t n_sites = sequences.front().size();

  out << '#';
  for (size_t i = 0; i < n_sequences; ++i)
    out << '\t' << node_names[i];
  out << '\n';

  for (size_t i = 0; i < n_sites; ++i) {
    out << i;
    for (size_t j = 0; j < n_sequences; ++j)
      out << '\t' << static_cast<bool>(sequences[j][i]);
    out << '\n';
  }
}

// struct watch_info {
//   double node_id;
//   double time_val;
//   double n_dom;
//   double mean_dom_size;
//   double dom_size_sd;
//   double fraction;
//   vector<size_t> patfreq;
// };


// std::ostream &
// operator<<(std::ostream &os, const watch_info &w) {
//   os << w.node_id << '\t'
//      << w.time_val << '\t'
//      << w.n_dom << '\t'
//      << w.mean_dom_size << '\t'
//      << w.dom_size_sd << '\t'
//      << w.fraction << '\t';
//   for (size_t ct = 0; ct < 8; ++ct)
//     os << w.patfreq[ct] << '\t';
//   return os;
// }


// static void
// update_watch_stats(const model_param &p,
//                    const PatSeq &patseq,
//                    const size_t node_id,
//                    const double watch,
//                    const double time,
//                    vector<watch_info> &w) {

//   w.push_back(watch_info());

//   vector<size_t> dom_size; // domain sizes
//   patseq.to_domain_sizes(dom_size);
//   const size_t n_domains = dom_size.size();

//   const double dom_tot =
//     std::accumulate(dom_size.begin(), dom_size.end(), 0.0);

//   const double mds = dom_tot/n_domains; // mean domain size

//   // standard deviation
//   const double sq_sum =
//     std::inner_product(dom_size.begin(), dom_size.end(), dom_size.begin(), 0.0);
//   const double stdev = std::sqrt(sq_sum/n_domains - mds*mds);

//   w.back().n_dom = n_domains;
//   w.back().mean_dom_size = mds;
//   w.back().dom_size_sd = stdev;
//   w.back().fraction = dom_tot/p.n_site;

//   patseq.get_all_context_freq(w.back().patfreq);

//   w.back().node_id = node_id;
//   w.back().time_val = floor(time/watch)*watch;
// }


////////////////////////////////////////////////////////////////////////////////
// SIMULATION
////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char **argv) {

  try {

    // static const size_t n_gibbs_iter = 500;

    string outfile;
    string pathfile;
    bool VERBOSE = false;
    bool EXTRA_VERBOSE = false;
    double watch = 0;
    string watchfile;
    size_t n_sites = 100;

    size_t rng_seed = std::numeric_limits<size_t>::max();

    OptionParser opt_parse(strip_path(argv[0]), "simulate methylome evolution",
                           "<params-file>");
    opt_parse.add_opt("output", 'o', "name of output file for methylomes"
                      "(default: stdout)", false, outfile);
    opt_parse.add_opt("n-sites", 'n', "length of sequence to simulate"
                      "(default: " + to_string(n_sites) + ")", false, n_sites);
    opt_parse.add_opt("paths", 'p', "name of output file for evolution paths"
                      "(default: stdout)", false, pathfile);
    opt_parse.add_opt("watch", 'w', "print summary statistics "
                      "at specified time interval (when -o)", false, watch);
    opt_parse.add_opt("seed", 's', "rng seed", false, rng_seed);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("extra-verbose", 'V', "print way more run info",
                      false, EXTRA_VERBOSE);
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
    const string param_file(leftover_args.back());
    if (!file_is_readable(param_file)) {
      cerr << "cannot read file: "<< param_file << endl
           << opt_parse.help_message() << endl;
      return  EXIT_SUCCESS;
    }
    ////////////////////////////////////////////////////////////////////////

    const bool keep_watch_info = (!pathfile.empty() && watch > 0);

    /* (1) INITIALIZING PARAMETERS */
    if (VERBOSE)
      cerr << "reading parameter file: " << param_file << endl;
    EpiEvoModel the_model;
    read_model(param_file, the_model);
    if (VERBOSE)
      cerr << the_model << endl;

    if (!pathfile.empty())
      write_pathfile_header(pathfile);

    /* standard mersenne_twister_engine seeded with rd()*/
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    if (VERBOSE)
      cerr << "rng seed: " << rng_seed << endl;
    std::mt19937 gen(rng_seed);

    /* (2) INITIALIZE THE ROOT SEQUENCE */
    if (VERBOSE)
      cerr << "[SIMULATING: " << the_model.node_names[0] << " (ROOT)]" << endl;
    vector<char> root_seq;
    the_model.sample_state_sequence_init(n_sites, gen, root_seq);

    StateSeq s(root_seq);
    if (VERBOSE) {
      cerr << "[SUMMARY:]" << endl
           << s.summary_string() << endl;
      vector<double> triplet_props;
      s.get_triplet_proportions(triplet_props);
      double total_rate = 0;
      for (size_t i = 0; i < triplet_props.size(); ++i)
        total_rate += the_model.triplet_rates[i]*triplet_props[i];
      // do we need to divide by the sequene length here?
      cerr << "mutations per site (at root): " << total_rate << endl;
    }

    /* EXTRACT INFO FROM THE PHYLOGENETIC TREE */
    const size_t n_nodes = the_model.t.get_size();

    vector<vector<char> > sequences(n_nodes, root_seq);
    // vector<watch_info> w; // only used if "watching" requested

    /* ITERATE OVER THE NODES IN THE TREE */
    for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
      const double curr_branch_len = the_model.branches[node_id];
      if (VERBOSE)
        cerr << "[SIMULATING: " << the_model.node_names[node_id] << " ("
             << curr_branch_len << ")]" << endl;

      TripletSampler ts(sequences[the_model.parent_ids[node_id]]);
      double time_value = 0;
      vector<Segment> the_path;

      /* SAMPLE CHANGES ALONG THE CURRENT BRANCH */
      while (time_value < curr_branch_len)
        sample_jump(the_model, curr_branch_len, gen, ts, the_path, time_value);

      ts.get_sequence(sequences[node_id]);

      if (!pathfile.empty())
        append_to_pathfile(pathfile, node_id, the_path);
    }

    if (!outfile.empty())
      write_output(outfile, the_model.node_names, sequences);

    // if (keep_watch_info) {
    //   watchfile = pathfile + ".stats";
    //   std::ofstream outstat(watchfile.c_str());
    //   outstat << "branch" << "\t" << "time" << "\t"
    //           << "n.domain" << "\t" << "mean.domain" << "\t"
    //           << "sd.domain" << "\t" << "fraction" << "\t"
    //           << "pattern000" << "\t" << "pattern001" << "\t"
    //           << "pattern010" << "\t" << "pattern011" << "\t"
    //           << "pattern100" << "\t" << "pattern101" << "\t"
    //           << "pattern110" << "\t" << "pattern111" << endl;
    //   copy(w.begin(), w.end(), std::ostream_iterator<watch_info>(outstat, "\n"));
    // }
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
