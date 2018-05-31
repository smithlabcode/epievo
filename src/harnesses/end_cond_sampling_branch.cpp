/* Copyright (C) 2018 University of Southern California
 *                    Jianghan Qu and Andrew D Smith
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
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_randist.h> /* chi-squared test */
#include <gsl/gsl_cdf.h>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "PhyloTreePreorder.hpp"
#include "Path.hpp"  /* related to Path */
#include "EpiEvoModel.hpp" /* model_param */
#include "TreeHelper.hpp"
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

static const size_t N_TRIPLETS = 8;
static const double MINWAIT = 1e-8;

////////////////////////////////////////////////////////////////////////////////
//////////   Record summary statistics                                //////////
////////////////////////////////////////////////////////////////////////////////

struct SummarySet {
  SummarySet(const Path l, const Path m, const Path r, const size_t n_bins);
  
  vector<size_t> num_jumps;
  vector<double> total_stay_time;
  vector<gsl_histogram *> h_time;
};

SummarySet::SummarySet(const Path l, const Path m, const Path r,
                       const size_t n_bins) {
  
  num_jumps.resize(N_TRIPLETS, 0);
  total_stay_time.resize(N_TRIPLETS, 0.0);
  
  vector<double> jumps(m.jumps.size() + 2, 0.0); // two extra entries
  copy(m.jumps.begin(), m.jumps.end(), jumps.begin() + 1); // start at jumps[1]
  jumps.back() = m.tot_time; // make sure the final entry is tot_time
  
  for(size_t k = 0; k < N_TRIPLETS; ++k) {
    gsl_histogram * h_time_k = gsl_histogram_alloc(n_bins);
    gsl_histogram_set_ranges_uniform(h_time_k, 0, m.tot_time+MINWAIT);
    h_time.push_back(h_time_k);
  }
  
  size_t context = triple2idx(l.init_state, m.init_state, r.init_state);
  for (size_t i = 1; i < jumps.size(); ++i) {
    const double t = 0.5*(jumps[i] + jumps[i-1]);
    context = triple2idx(l.state_at_time(t), m.state_at_time(t),
                         r.state_at_time(t));
    ++num_jumps[context];
    total_stay_time[context] += (jumps[i] - jumps[i-1]);
    gsl_histogram_increment(h_time[context], jumps[i] - jumps[i-1]);
    
  }
  --num_jumps[context]; // last break point is not a jump
}

////////////////////////////////////////////////////////////////////////////////
//////////   Report and Test summaries                                //////////
////////////////////////////////////////////////////////////////////////////////
static double
evaluate_fit(const vector<double> &reference,
             const vector<double> &to_evaluate) {
  
  assert(reference.size() == to_evaluate.size());
  
  double chi_squared_stat = 0.0;
  for (size_t i = 0; i < reference.size(); ++i) {
    chi_squared_stat +=
    (to_evaluate[i] - reference[i])*
    (to_evaluate[i] - reference[i])/reference[i];
  }
  
  const double degrees_of_freedom = reference.size() - 1;
  return gsl_cdf_chisq_P(chi_squared_stat, degrees_of_freedom);
}


struct SummaryStatsFreq {
  SummaryStatsFreq(const vector<SummarySet> &summary);
  string print_time(const size_t triplet_idx) const;
  string print_jumps(const size_t triplet_idx) const;
  
  size_t num_samples;
  vector<vector<double> > time_freq; // proportion, not count
  vector<vector<double> > jumps_freq; // proportion, not count
};

SummaryStatsFreq::SummaryStatsFreq(const vector<SummarySet> &summary) {
  num_samples = summary.size();
  
  if (num_samples > 0) {
    for(size_t k = 0; k < N_TRIPLETS; ++k) {
      vector<double> time_freq_k;
      vector<double> jumps_freq_k;
      
      // initialize time_freq
      const size_t nbins_time = summary.back().h_time[k]->n;
      time_freq_k.resize(nbins_time, 0);
      
      // initialize jumps_freq
      size_t max_jumps = 0;
      for (size_t i = 0; i < num_samples; i++)
        if (summary[i].num_jumps[k] > max_jumps)
          max_jumps = summary[i].num_jumps[k];
      jumps_freq_k.resize(max_jumps + 1, 0);
      
      // merge summaries
      for (size_t i = 0; i < num_samples; i++) {
        for (size_t j = 0; j < nbins_time; j++) {
          time_freq_k[j] += summary[i].h_time[k]->bin[j];
        }
        jumps_freq_k[summary[i].num_jumps[k]]++;
      }
      
      // normalize to frequency
      const double sum_time = accumulate(time_freq_k.begin(), time_freq_k.end(),
                                         0.0);
      const double sum_jumps = accumulate(jumps_freq_k.begin(),
                                          jumps_freq_k.end(), 0.0);
      
      for (size_t i = 0; i < time_freq_k.size(); i++)
        time_freq_k[i] /= sum_time;
      
      for (size_t i = 0; i < jumps_freq_k.size(); i++)
        jumps_freq_k[i] /= sum_jumps;
    
      time_freq.push_back(time_freq_k);
      jumps_freq.push_back(jumps_freq_k);
    }
  }
}

string
SummaryStatsFreq::print_time(const size_t triplet_idx) const {
  string str;
  if (time_freq[triplet_idx].size() > 0)
    str = std::to_string(time_freq[triplet_idx][0]);
  for(size_t i = 1; i < time_freq[triplet_idx].size(); i++) {
    str = str + "," + std::to_string(time_freq[triplet_idx][i]);
  }
  return str;
}


string
SummaryStatsFreq::print_jumps(const size_t triplet_idx) const {
  string str;
  if (jumps_freq[triplet_idx].size() > 0)
    str = std::to_string(jumps_freq[triplet_idx][0]);
  for(size_t i = 1; i < jumps_freq[triplet_idx].size(); i++) {
    str = str + "," + std::to_string(jumps_freq[triplet_idx][i]);
  }
  return str;
}

static void
test_summary(SummaryStatsFreq &a, SummaryStatsFreq &b,
             double &pval_time, double &pval_jumps,
             const size_t triplet_idx) {
  // rule out zero expected values
  vector<double> time_exp, time_obs, jump_exp, jump_obs;
  for(size_t i = 0; i < std::min(a.time_freq[triplet_idx].size(),
                                 b.time_freq[triplet_idx].size()); i++) {
    if (a.time_freq[triplet_idx][i] > 0) {
      time_exp.push_back(a.time_freq[triplet_idx][i]);
      time_obs.push_back(b.time_freq[triplet_idx][i]);
    }
  }
  pval_time = evaluate_fit(time_exp, time_obs);
  
  for(size_t i = 0; i < std::min(a.jumps_freq[triplet_idx].size(),
                                 b.jumps_freq[triplet_idx].size()); i++) {
    if (a.jumps_freq[triplet_idx][i] > 0) {
      jump_exp.push_back(a.jumps_freq[triplet_idx][i]);
      jump_obs.push_back(b.jumps_freq[triplet_idx][i]);
    }
  }
  pval_jumps = evaluate_fit(jump_exp, jump_obs);
}



////////////////////////////////////////////////////////////////////////////////
//////////   copied/modified from SingleSampler.cpp                   //////////
////////////////////////////////////////////////////////////////////////////////
// collect rates and interval lengths
static void
collect_segment_info(const vector<double> &triplet_rates,
                     const Path &l, const Path &r,
                     vector<SegmentInfo> &seg_info) {
  Environment env(l, r);
  const size_t n_intervals = env.left.size();
  seg_info = vector<SegmentInfo>(n_intervals);
  
  for (size_t i = 0; i < n_intervals; ++i) {
    const size_t pattern0 = triple2idx(env.left[i], false, env.right[i]);
    const size_t pattern1 = triple2idx(env.left[i], true, env.right[i]);
    seg_info[i] = SegmentInfo(triplet_rates[pattern0], triplet_rates[pattern1],
                              env.breaks[i] - (i == 0 ? 0.0 : env.breaks[i-1]));
    assert(seg_info[i].len > 0.0);
  }
}

static double
root_post_prob0(const size_t site_id, const vector<Path> &the_paths,
                const vector<vector<double> > &horiz_tr_prob,
                const vector<double> &q) {
  
  const size_t left_st = the_paths[site_id - 1].init_state;
  const size_t right_st = the_paths[site_id + 1].init_state;
  
  const double p0 = (horiz_tr_prob[left_st][0]*horiz_tr_prob[0][right_st])*q[0];
  const double p1 = (horiz_tr_prob[left_st][1]*horiz_tr_prob[1][right_st])*q[1];
  
  return p0/(p0 + p1);
}

////////////////////////////////////////////////////////////////////////////////
//////////   Downward_sampling_branch Forward sampling                //////////
////////////////////////////////////////////////////////////////////////////////
static void
forward_sample(const vector<double> &rates,
               const bool is, const bool es,
               const double tot_time, std::mt19937 &gen,
               vector<double> &jump_times, const size_t max_iteration) {
  
  bool reach_target = false;
  vector<double> fs_jump_times;
  size_t num_sampled = 0;
  
  while (!reach_target && num_sampled < max_iteration) {
    fs_jump_times.clear();
    double time_value = 0;
    size_t curr_state = is;
    
    // one forward sampling run
    while (time_value < tot_time) {
      const double holding_rate = rates[curr_state];
      std::exponential_distribution<double> exp_distr(holding_rate);
      const double holding_time =
      std::max(exp_distr(gen), std::numeric_limits<double>::min());
      
      time_value += holding_time;
      
      if (time_value < tot_time) {
        fs_jump_times.push_back(time_value);
        curr_state = complement_state(curr_state);
      }
    }
    reach_target = curr_state == es;
    ++num_sampled;
  }
  
  if (reach_target)
    for (size_t i = 0; i < fs_jump_times.size(); ++i) {
      jump_times.push_back(fs_jump_times[i]);
    }
}


// using Forward sampling + rejection
static void
downward_sampling_branch_fs(const size_t site, const vector<Path> &paths,
                            const vector<SegmentInfo> &seg_info,
                            const FelsHelper &fh,
                            std::mt19937 &gen, Path &new_path,
                            const size_t max_iterations) {
  
  const size_t n_intervals = seg_info.size();
  
  size_t par_state = new_path.init_state;
  double time_passed = 0;
  
  for (size_t m = 0; m < n_intervals - 1; ++m) {
    // compute conditional posterior probability
    vector<vector<double> > P; // transition prob matrix
    continuous_time_trans_prob_mat(seg_info[m].rate0, seg_info[m].rate1,
                                   seg_info[m].len, P);
    
    const double p0 = P[par_state][0] / fh.p[m][par_state] *
                      ((m == n_intervals - 1) ? fh.q[0] : fh.p[m + 1][0]);
    
    // generate random state at break point
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    bool new_state = (unif(gen) > p0);
    
    // generate path
    vector<double> jump_times;
    const vector<double> rates {seg_info[m].rate0, seg_info[m].rate1};
    forward_sample(rates, par_state, new_state, seg_info[m].len, gen,
                   jump_times, max_iterations);
    
    // append jump_times to new_path
    for (size_t i = 0; i < jump_times.size(); ++i) {
      new_path.jumps.push_back(time_passed + jump_times[i]) ;
    }
    // prepare for next interval
    time_passed += seg_info[m].len;
    par_state = new_state;
  }
  
  // we only test root-leaf branch.
  const bool leaf_state = paths[site].end_state();
  // generate path
  vector<double> jump_times;
  const size_t m = n_intervals - 1;
  // prepare helper values
  const vector<double> rates {seg_info[m].rate0, seg_info[m].rate1};
  forward_sample(rates, par_state, leaf_state, seg_info[m].len, gen,
                 jump_times, max_iterations);
  // append jump_times to new_path
  for (size_t i = 0; i < jump_times.size(); ++i) {
    new_path.jumps.push_back(time_passed + jump_times[i]) ;
  }
  
  assert(new_path.end_state() == leaf_state);
}


////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char **argv) {
  try {

    bool VERBOSE = false;
    string outfile;
    string outstatefile;

    size_t max_iterations = 1000000;
    size_t n_paths_to_sample = 1000;
    size_t n_hist_time_bins = 5;
    size_t test_site = 2; // we test 5-site path
    size_t test_branch = 1; // branch to test in the tree
    
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
    opt_parse.add_opt("paths", 'n', "number of paths to sample",
                      false, n_paths_to_sample);
    opt_parse.add_opt("bins", 'b', "number of bins of holding time histogram",
                      false, n_hist_time_bins);
    opt_parse.add_opt("outfile", 'o', "outfile (sampling statistics)",
                      true, outfile);
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
      cerr << "[READING PARAMETER FILE: " << param_file << ", "
      << tree_file << "]" << endl;
    
    EpiEvoModel the_model;
    read_model(param_file, the_model);
    if (VERBOSE)
      cerr << the_model << endl;
    
    if (VERBOSE)
      cerr << "[READING TREE FILE: " << tree_file << "]" << endl;
    PhyloTreePreorder the_tree; // tree topology and branch lengths
    std::ifstream tree_in(tree_file.c_str());
    if (!tree_in || !(tree_in >> the_tree))
      throw std::runtime_error("bad tree file: " + tree_file);
    const size_t n_nodes = the_tree.get_size();
    TreeHelper th(the_tree);

    if (VERBOSE)
      cerr << "[READING PATHS: " << pathsfile << "]" << endl;
    vector<vector<Path> > paths; // along multiple branches
    vector<string> node_names;
    read_paths(pathsfile, node_names, paths);

    //const size_t n_nodes = node_names.size();
    /* below: 1st element of paths empty at root; use last */
    //const size_t n_sites = paths.back().size();

    if (VERBOSE)
      cerr << "TEST BRANCH: " << test_branch << endl
           << "TEST SITE: " << test_site << endl
           << "PATHS TO SIMULATE: " << n_paths_to_sample << endl;
    // standard mersenne_twister_engine seeded with rd()
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    if (VERBOSE)
      cerr << "rng seed: " << rng_seed << endl;

    std::mt19937 gen(rng_seed);

    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    cerr << "----- TEST upward downward sampling BELOW ---------" << endl;
  
    vector<vector<SegmentInfo> > seg_info(n_nodes);
    for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
      collect_segment_info(the_model.triplet_rates,
                           paths[node_id][test_site - 1],
                           paths[node_id][test_site + 1], seg_info[node_id]);
    }
    
    vector<SummarySet> fs_summary, ds_summary;
    size_t progress = 0;
    
    // (1) Upward pruning
    vector<FelsHelper> fh;
    pruning(th, test_site, paths, seg_info, fh);
    
    while (fs_summary.size() < n_paths_to_sample) {
      if (VERBOSE && fs_summary.size() * 10 % n_paths_to_sample == 0) {
        cerr << "FINISHED: " << progress << '%' << endl;
        progress += 10;
      }
      
      // (2) Downward sampling: Direct sampling
      vector<Path> new_path_ds_all;
      downward_sampling(th, test_site, paths, the_model.init_T, seg_info,
                        fh, gen, new_path_ds_all);
      
      SummarySet current_summary_ds(paths[test_branch][test_site-1],
                                    new_path_ds_all[test_branch],
                                    paths[test_branch][test_site+1],
                                    n_hist_time_bins);
      ds_summary.push_back(current_summary_ds);
      
      // (3) Dowanward sampling: forward sampling
      Path new_path_fs;
      
      // sample new root state
      const double root_p0 = root_post_prob0(test_site, paths[0],
                                             the_model.init_T, fh[0].q);

      std::uniform_real_distribution<double> unif(0.0, 1.0);
      bool new_root_state = (unif(gen) > root_p0);
      new_path_fs.init_state = new_root_state;
      new_path_fs.tot_time = paths[test_branch][test_site - 1].tot_time;
      
      // preorder traversal of the tree
      downward_sampling_branch_fs(test_site, paths[test_branch],
                                  seg_info[test_branch], fh[test_branch], gen,
                                  new_path_fs, max_iterations);
      SummarySet current_summary_fs(paths[test_branch][test_site-1],
                                    new_path_fs,
                                    paths[test_branch][test_site+1],
                                    n_hist_time_bins);
      fs_summary.push_back(current_summary_fs);
    }
    cerr << "FINISHED: " << progress << '%' << endl;
    
    // get the summaries of the summaries
    SummaryStatsFreq FS_report(fs_summary);
    SummaryStatsFreq DS_report(ds_summary);
    
    // write output
    std::ofstream out(outfile.c_str());
    out << "TEST_SITE_BRANCH" << '\t' << test_site << '\t'
                              << test_branch << endl;
    out << "N_PATHS"    << '\t' << n_paths_to_sample << endl;
    out << "IN_STATE" << '\t'
        << "PVAL_JUMPS" << '\t' << "PVAL_TIME" << '\t'
        << "HIST_JUMPS_FS" << '\t' << "HIST_JUMPS_DS" << '\t'
        << "HIST_TIME_FS" << '\t' << "HIST_TIME_DS" << endl;
    
    for(size_t k = 0; k < N_TRIPLETS; ++k) {
      double pval_time, pval_jumps;
      test_summary(FS_report, DS_report, pval_time, pval_jumps, k);
      out << std::bitset<3>(k).to_string() << '\t'
          << pval_jumps << '\t' << pval_time << '\t'
          << FS_report.print_jumps(k) << '\t'
          << DS_report.print_jumps(k) << '\t'
          << FS_report.print_time(k) << '\t'
          << DS_report.print_time(k) << endl;
    }
    

  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
