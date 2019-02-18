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
#include <cassert>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <random>
#include <functional>
#include <iomanip>

#include "IndepSite.hpp"
#include "PhyloTreePreorder.hpp"
#include "Path.hpp"

using std::vector;
using std::string;


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
               const vector<vector<Path> > &paths, const vector<double> &rates,
               vector<FelsHelper> &fh) {
  
  // avoid recursion by iterating backwards ensuring all child nodes
  // are processed before their parents
  fh.resize(th.n_nodes);
  for (size_t i = th.n_nodes; i > 0; --i) {
    const size_t node_id = i - 1;
    
    vector<double> q(2, 1.0);
    if (th.is_leaf(node_id)) {
      const bool leaf_state = paths[site_id].end_state();
      q[0] = (leaf_state == false) ? 1.0 : 0.0;
      q[1] = (leaf_state == true)  ? 1.0 : 0.0;
    }
    else {
      for (ChildSet c(th.subtree_sizes, node_id); c.good(); ++c) {
        q[0] *= fh[*c].p.front()[0];
        q[1] *= fh[*c].p.front()[1];
      }
    }
    fh[node_id].q.swap(q); // assign computed q to the fh
    
    /* now calculate p if we are not at root */
    if (!th.is_root(node_id)) {
      vector<vector<double> > P;
      continuous_time_trans_prob_mat(rates[0], rates[1],
                                     paths[node_id][site_id].tot_time, P);
      // p <- P*q
      const double p = { P[0][0]*q[0] + P[0][1]*q[1],
        P[1][0]*q[0] + P[1][1]*q[1] };
      fh[node_id].p.swap(p); // assign computed p to the fh
    }
  }
}


static double
root_post_prob0(const size_t site_id, const vector<double> init_pi,
                const vector<double> &q) {
  
  const double p0 = init_pi[0]*q[0];
  const double p1 = init_pi[1]*q[1];
  
  return p0/(p0 + p1);
}


/* Joint posterior distribution p(u, v | X(u)) */
static void
joint_post(const double p0_u, const FelsHelper &fh,
           const vector<vector<double> > &PT,
           vector<vector<double> > &p_joint) {
  
  p_joint[0][0] = PT[0][0] * fh.q[0] * p0_u / fh.p[0];
  p_joint[0][1] = PT[0][1] * fh.q[1] * p0_u / fh.p[0];
  p_joint[1][0] = PT[1][0] * fh.q[0] * (1 - p0_u) / fh.p[1];
  p_joint[1][1] = PT[1][1] * fh.q[1] * (1 - p0_u) / fh.p[1];
  
  const double Z = p_joint[0][0] + p_joint[0][1] + p_joint[1][0] + p_joint[1][1];
  assert(Z > 0);
  
  for (size_t u = 0; u < 2; u++)
    for (size_t v = 0; v < 2; v++)
      p_joint[u][v] /= Z;
}


static void
weighted_J_D_branch(const vector<double> &rates, const double T,
                    const FelsHelper fh,
                    const double p0_parent, double &p0_child,
                    vector<double> &J, vector<double> &D) {

  vector<vector<double> > P;
  continuous_time_trans_prob_mat(rates[0], rates[1], T, P);
  
  // get joint distribution of end states
  vector<vector<double> > p_joint(2, vector<double> (2, 0.0));
  joint_post(p0_parent, P, p_joint);
  p0_child = p_joint[0][0] + p_joint[1][0];
  
  // compute weighted expectation of sufficient stats
  vector< vector<double>> J0(2, vector<double> (2, 0.0));
  vector< vector<double>> J1(2, vector<double> (2, 0.0));
  vector< vector<double>> D0(2, vector<double> (2, 0.0));
  vector< vector<double>> D1(2, vector<double> (2, 0.0));
  expectation_J(rates, paths[node_id][site_id].tot_time, J0, J1);
  expectation_D(rates, paths[node_id][site_id].tot_time, D0, D1);
  
  for (size_t u = 0; u < 2; u++)
    for (size_t v = 0; v < 2; v++) {
      J[0] += p_joint[u][v]*J0[u][v];
      J[1] += p_joint[u][v]*J1[u][v];
      D[0] += p_joint[u][v]*D0[u][v];
      D[1] += p_joint[u][v]*D1[u][v];
  }
}



/* Iterates over all nodes in pre-order, and for each branch it computes
 * conditional means of sufficient statistics.
 */
static void
pruning_downward(const TreeHelper &th, const size_t site_id,
                 const vector<vector<Path> > &paths,
                 const vector<double> init_pi, const vector<double> &rates,
                 const vector<FelsHelper> &fh,
                 vector<vector<double> > &J,
                 vector<vector<double> > &D) {
  
  vector<double> p0_margin(node_id);
  
  // marginalize at root node
  p0_margin[0] = root_post_prob0(site_id, init_pi, fh[0].q);
  
  for (size_t node_id = 1; node_id < th.n_nodes; ++node_id)
    weighted_J_D_branch(rates, T, fh[node_id],
                        p0_margin[th.parent_ids[node_id]],
                        p0_margin[node_id], J[node_id], D[node_id]);
}


/* Iterates over all nodes in pre-order, and for each branch it samples
 * a new path.
 */
static void
sampling_downward(const TreeHelper &th, const size_t site_id,
                  const vector<vector<Path> > &paths,
                  const vector<double> init_pi, const vector<double> &rates,
                  const vector<FelsHelper> &fh,
                  std::mt19937 &gen,
                  vector<Path> &sampled_path) {
  
  vector<double> p0_margin(node_id);
  std::uniform_real_distribution<double> unif(0.0, 1.0);

  // marginalize at root node
  p0_margin[0] = root_post_prob0(site_id, init_pi, fh[0].q);
  sampled_path = vector<Path>(th.n_nodes);
  sampled_path.front() = Path(unif(gen) > p0_margin[0], 0.0);

  for (size_t node_id = 1; node_id < th.n_nodes; ++node_id) {
    const size_t start_state = sampled_path[th.parent_ids[node_id]].end_state();
    const double T = paths[node_id][site_id].tot_time;

    sampled_path[node_id].init_state = start_state;
    sampled_path[node_id].tot_time = T;
    
    vector<vector<double> > P;
    continuous_time_trans_prob_mat(rates[0], rates[1], T, P);
    const double p0 =
    P[start_state][0] * fh[node_id].q[0] / fh[node_id].p[start_state];
    const size_t sampled_state = (unif(gen) > p0);
    
    const CTMarkovModel ctmm(seg_info[i].rate0, seg_info[i].rate1);
    assert(end_cond_sample_forward_rejection(10000000, ctmm,
                                             start_state, sampled_state,
                                             T, gen, sampled_path.jumps, 0));
  }
}

////////////////////////////////////////////////////////////////////////////////
///////////////      SUFF-STATS CONDITIONAL MEANS     //////////////////////////
////////////////////////////////////////////////////////////////////////////////

static double
expectation_J(const vector<double> rates, const double T,
              vector< vector<double>> &J0, vector< vector<double>> &J1) {
  const double r0 = rates[0];
  const double r1 = rates[1];
  const double s = r0 + r1;
  const double p = r0 * r1;
  const double d = r1 - r0;
  const double e = exp(- s * T);
  
  const double C1 = d * (1 - e) / s;
  J0[0][0] = p * ( T * (r1 - r0*e) - C1 ) / (s * (r1 + r0*e));
  J1[0][0] = J0[0][0];
  
  J0[1][1] = p * ( T * (r0 - r1*e) + C1 ) / (s * (r0 + r1*e));
  J1[1][1] = J0[1][1];
  
  const double C2 = p * T * (1 + e) / (s * (1 - e));
  const double C3 = (r0 * r0 + r1 * r1) / (s * s);
  const double C4 = (2 * p) / (s * s);
  
  J0[0][1] = C2 + C3;
  J1[0][1] = C2 - C4;
  J0[1][0] = J1[0][1];
  J1[1][0] = J0[0][1];
}


static double
expectation_D(const vector<double> rates, const double T,
              vector< vector<double>> &D0, vector< vector<double>> &D1) {
  const double r0 = rates[0];
  const double r1 = rates[1];
  const double r00 = r0 * r0;
  const double r11 = r1 * r1;
  const double s = r0 + r1;
  const double p = r0 * r1;
  const double e = exp(- s * T);
  
  const double C1 = 2 * p * (1 - e) / s;
  const double C2 = p * T * (1 + e);
  D0[0][0] = ((r11 + r00 * e) * T + C1) / (s * (r1 + r0 * e));
  D1[0][0] = (C2 * T - C1) / (s * (r1 + r0 * e));
  
  D0[1][1] = (C2 * T - C1) / (s * (r0 + r1 * e));
  D1[1][1] = ((r00 + r11 * e) * T + C1) / (s * (r0 + r1 * e));
  
  const double C3 = (p - r00) * (1 - e) / s;
  D0[0][1] = ((p - r00 * e) * T - C3) / (s * (r0 - r0 * e));
  D1[0][1] = ((r00 - p * e) * T + C3) / (s * (r0 - r0 * e));
  
  const double C4 = (p - r11) * (1 - e) / s;
  D0[0][1] = ((r11 - p * e) * T + C4) / (s * (r1 - r1 * e));
  D1[0][1] = ((p - r11 * e) * T - C4) / (s * (r1 - r1 * e));
}


////////////////////////////////////////////////////////////////////////////////
///////////////          site-independent E-step       /////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* Obtain conditional means of sufficient statistics */
void
expectation_sufficient_statistics(const vector<double> rates,
                                  const vector<double> init_pi,
                                  const TreeHelper &th,
                                  const vector<vector<Path> > &paths,
                                  vector<vector<double> > &J,
                                  vector<vector<double> > &D) {
  
  assert(paths.size() > 1);
  
  /* sufficient statistics: J - jump counts, D - spent time */
  J.resize(th.n_nodes);
  D.resize(th.n_nodes);
  for (size_t node_id = 1; node_id < th.n_nodes; ++node_id) {
    J[node_id].resize(2);
    J[node_id].resize(2);
  }
  
  
  for (size_t site_id = 0; site_id < paths[1].size(); site_id++) {
    vector<FelsHelper> fh(th.n_nodes);
    pruning_upward(th, site_id, paths, rates, fh);
    pruning_downward(th, site_id, paths, init_pi, rates, fh, J, D);
  }
}


/* Sample new evolution history */
void
sample_paths(const vector<double> rates, const vector<double> init_pi,
             const TreeHelper &th,
             const vector<vector<Path> > &paths,
             std::mt19937 &gen, vector<vector<Path> > &sampled_paths) {
  
  assert(paths.size() > 1);
  sampled_paths.resize(paths[1].size());
  
  for (size_t site_id = 0; site_id < paths[1].size(); site_id++) {
    vector<FelsHelper> fh(th.n_nodes);
    pruning_upward(th, site_id, paths, rates, fh);
    sampling_downward(th, site_id, paths, init_pi, rates, fh,
                      gen, sampled_paths[site_id]);
  }
}

////////////////////////////////////////////////////////////////////////////////
///////////////          site-independent M-step       /////////////////////////
////////////////////////////////////////////////////////////////////////////////

void
compute_estimates_rates_and_branches(const vector<vector<double> > &J,
                                     const vector<vector<double> > &D,
                                     vector<double> &rates, TreeHelper &th) {
  
  vector<double> updated_rates;
  vector<double> updated_branches;
  
  vector<double> J_sum(2, 0.0);
  vector<double> D_sum(2, 0.0);

  for (size_t b = 1; b < branches.size(); ++b) {
    J_sum[0] += J[b][0];
    J_sum[1] += J[b][1];

    D_sum[0] += D[b][0];
    D_sum[1] += D[b][1];
  }
  
  // update rates
  assert(D_sum[0] > 0 && D_sum[1] > 0);
  rates[0] = J_sum[0] / D_sum[0];
  rates[1] = J_sum[1] / D_sum[1];

  // update branch lengths
  for (size_t b = 1; b < branches.size(); ++b) {
    th.branches[b] = (J[b][0]*rates[0] + J[b][1]*rates[1]) / (D[b][0] + D[b][1]);
  }
  
  /* we don't do any scaling here */
}


void
estimate_root_distribution(const vector<vector<Path> > &paths,
                           vector<double> &init_pi) {
  
  assert(paths.size() > 1 && !paths[1].empty());
  
  init_pi.resize(2);

  double N1 = 0.0;
  for(size_t site_id = 0; site_id < paths[1].size(); site_id++)
    N1 += paths[1][site_id].init_state;

  init_pi[1] = N1 / paths[1].size();
  init_pi[0] = 1 - init_pi[1];
}
