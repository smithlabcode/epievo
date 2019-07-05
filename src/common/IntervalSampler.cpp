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
#include <cassert>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <random>
#include <functional>
#include <iomanip>

#include "IntervalSampler.hpp"
#include "PhyloTreePreorder.hpp"
#include "Path.hpp"
#include "StateSeq.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::pair;
using std::make_pair;

using std::bind;
using std::placeholders::_1;
using std::multiplies;
using std::runtime_error;


/* collect rates and interval lengths */
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


////////////////////////////////////////////////////////////////////////////////
//////////////////////              PROPOSAL           /////////////////////////
////////////////////////////////////////////////////////////////////////////////

static void
propose_h_interval(const vector<vector<Path> > &paths,
                   const size_t site_id, const size_t node_id,
                   Path &proposed_path, const vector<SegmentInfo> &seg_info,
                   std::mt19937 &gen, const double time_to_update,
                   double &orig_proposal, double &update_proposal) {
  
  // scan and copy path
  proposed_path = Path(paths[node_id][site_id].init_state,
                       paths[node_id][site_id].tot_time);
  size_t start_state = paths[node_id][site_id].init_state;

  size_t start_jump = 0;
  double start_time = 0.0;
  double end_time = seg_info[0].len;
  size_t seg_id = 0;
  
  // copy the jumps above interested segment and determine start state
  while (end_time < time_to_update && seg_id < seg_info.size() - 1) {
    // start_jump is the first jump in interested interval
    while (start_jump < paths[node_id][site_id].jumps.size() &&
           paths[node_id][site_id].jumps[start_jump] < end_time) {
      proposed_path.jumps.push_back(paths[node_id][site_id].jumps[start_jump]);
      start_state = complement_state(start_state);
      start_jump++;
    }
    
    start_time = end_time;
    end_time += seg_info[++seg_id].len;
  }
  
  // determine the end state and the first jump after interested interval
  size_t end_jump = start_jump;
  size_t end_state = start_state;
  while (end_jump < paths[node_id][site_id].jumps.size() &&
         paths[node_id][site_id].jumps[end_jump] < end_time)
    end_jump++;
  if ((end_jump - start_jump) % 2 == 1)
    end_state = complement_state(end_state);
  
  // end-conditioned sampling
  const TwoStateCTMarkovModel ctmm(seg_info[seg_id].rate0,
                                   seg_info[seg_id].rate1);
  end_cond_sample_Poisson(ctmm, start_state, end_state,
                          seg_info[seg_id].len, gen, proposed_path.jumps,
                          start_time);
  size_t end_jump_update = proposed_path.jumps.size();
  
  // calculate proposal probabilities
  orig_proposal =
  end_cond_sample_Poisson_prob(ctmm, paths[node_id][site_id].jumps,
                               start_state, end_state,
                               start_time, end_time, start_jump, end_jump);
  update_proposal =
  end_cond_sample_Poisson_prob(ctmm, proposed_path.jumps, start_state, end_state,
                               start_time, end_time,
                               start_jump, end_jump_update);
  
  // copy rest jumps to proposed path
  while (end_jump < paths[node_id][site_id].jumps.size())
    proposed_path.jumps.push_back(paths[node_id][site_id].jumps[end_jump++]);
  
  assert(proposed_path.init_state == paths[node_id][site_id].init_state);
  assert(proposed_path.end_state() == paths[node_id][site_id].end_state());
}



////////////////////////////////////////////////////////////////////////////////
///////////////        Compute acceptance rate         /////////////////////////
////////////////////////////////////////////////////////////////////////////////
static double
log_likelihood(const vector<double> &rates,
               const vector<double> &J, const vector<double> &D) {
  static const size_t n_triples = 8;
  double r = 0.0;
  for (size_t i = 0; i < n_triples; ++i) {
    r += J[i]*log(rates[i]) - D[i]*rates[i];
  }
  return r;
}


/* compute acceptance rate */
double
log_accept_rate(const EpiEvoModel &mod,
                const size_t site_id, const size_t node_id,
                const vector<vector<Path> > &paths,
                const Path &proposed_path,
                const double orig_proposal, const double update_proposal) {
  
  // calculate likelihood
  static const size_t n_triples = 8;
  
  double llr = orig_proposal - update_proposal;

  //cerr << "orig_jumps: " << paths[node_id][site_id].jumps.size()
  //<< ", update_jumps: " << proposed_path.jumps.size() << endl;
  //cerr << "orig_prop: " << orig_proposal
  //<< ", update_prop: " << update_proposal << endl;
  //cerr << "(orig_proposal - update_proposal): " << llr << endl;
  
  vector<double> D_orig(n_triples), J_orig(n_triples);
  fill_n(begin(D_orig), n_triples, 0.0);
  fill_n(begin(J_orig), n_triples, 0.0);

  vector<double> D_prop(n_triples), J_prop(n_triples);
  fill_n(begin(D_prop), n_triples, 0.0);
  fill_n(begin(J_prop), n_triples, 0.0);
  
  if (site_id > 1) {
    add_sufficient_statistics(paths[node_id][site_id - 2],
                              paths[node_id][site_id - 1],
                              paths[node_id][site_id], J_orig, D_orig);
    add_sufficient_statistics(paths[node_id][site_id - 2],
                              paths[node_id][site_id - 1],
                              proposed_path, J_prop, D_prop);
  }

  add_sufficient_statistics(paths[node_id][site_id - 1],
                            paths[node_id][site_id],
                            paths[node_id][site_id + 1], J_orig, D_orig);
  add_sufficient_statistics(paths[node_id][site_id - 1],
                            proposed_path,
                            paths[node_id][site_id + 1], J_prop, D_prop);
  
  
  if (site_id < paths[node_id].size() - 2) {
    add_sufficient_statistics(paths[node_id][site_id],
                              paths[node_id][site_id + 1],
                              paths[node_id][site_id + 2], J_orig, D_orig);
    add_sufficient_statistics(proposed_path,
                              paths[node_id][site_id + 1],
                              paths[node_id][site_id + 2], J_prop, D_prop);
  }

  llr += (log_likelihood(mod.triplet_rates, J_prop, D_prop) -
          log_likelihood(mod.triplet_rates, J_orig, D_orig));

  //cerr << "llr: " << llr << endl;
  //cerr << endl << endl;
  return llr;
}



// Currently this function only works for pairwise inference.
bool
Metropolis_Hastings_interval(const EpiEvoModel &the_model, const TreeHelper &th,
                             const size_t site_id, vector<vector<Path> > &paths,
                             std::mt19937 &gen) {
  
  // randomly pick a node
  std::discrete_distribution<size_t> node_sampler(th.branches.begin(),
                                                  th.branches.end());
  const size_t node_to_update = node_sampler(gen);
  vector<SegmentInfo> seg_info;
  collect_segment_info(the_model.triplet_rates,
                       paths[node_to_update][site_id - 1],
                       paths[node_to_update][site_id + 1], seg_info);
  
  // randomly pick a time
  std::uniform_real_distribution<double> time_sampler(0.0,
                                                      paths[node_to_update][site_id].tot_time);
  const double time_to_update = time_sampler(gen);

  // propose segment
  double orig_proposal, update_proposal;
  Path proposed_path;
  propose_h_interval(paths, site_id, node_to_update, proposed_path, seg_info,
                     gen, time_to_update, orig_proposal, update_proposal);
  

  // acceptance rate
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  const double u = unif(gen) ;
  const double log_acc_rate =
  log_accept_rate(the_model, site_id, node_to_update, paths, proposed_path,
                  orig_proposal, update_proposal);
  
  bool accepted = false;
  if (log_acc_rate >= 0 || u < exp(log_acc_rate))
    accepted = true;
  
  // if accepted, replace old path with proposed one.
  if (accepted) 
    paths[node_to_update][site_id] = proposed_path;

  return accepted;
}
