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
#include "EpiEvoModel.hpp"
#include "ContinuousTimeMarkovModel.hpp"
#include "EndCondSampling.hpp"


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
upward_process(const TreeHelper &th, const size_t site_id,
               const vector<vector<Path> > &paths, const vector<double> &rates,
               vector<FelsHelper> &fh) {

  // avoid recursion by iterating backwards ensuring all child nodes
  // are processed before their parents
  fh.resize(th.n_nodes);
  for (size_t i = th.n_nodes; i > 0; --i) {
    const size_t node_id = i - 1;

    vector<double> q(2, 1.0);
    if (th.is_leaf(node_id)) {
      const bool leaf_state = paths[site_id][node_id].end_state();
      q[0] = (leaf_state == false) ? 1.0 : 0.0;
      q[1] = (leaf_state == true)  ? 1.0 : 0.0;
    }
    else {
      for (ChildSet c(th.subtree_sizes, node_id); c.good(); ++c) {
        q[0] *= fh[*c].p[0];
        q[1] *= fh[*c].p[1];
      }
    }
    fh[node_id].q.swap(q); // assign computed q to the fh

    /* now calculate p if we are not at root */
    if (!th.is_root(node_id)) {
      two_by_two P;
      continuous_time_trans_prob_mat(rates[0], rates[1],
                                     paths[site_id][node_id].tot_time, P);
      // p <- P*q
      vector<double> p = { P[0][0]*fh[node_id].q[0] + P[0][1]*fh[node_id].q[1],
        P[1][0]*fh[node_id].q[0] + P[1][1]*fh[node_id].q[1] };
      fh[node_id].p.swap(p); // assign computed p to the fh
    }
  }
}


static double
root_post_prob0(const vector<double> &init_pi, const vector<double> &q) {

  const double p0 = init_pi[0]*q[0];
  const double p1 = init_pi[1]*q[1];

  return p0/(p0 + p1);
}


/* Joint posterior distribution p(u, v | X(u)) */
static void
joint_post(const double p0_u, const FelsHelper &fh,
           const two_by_two &PT,
           two_by_two &p_joint) {

  p_joint(0, 0) = PT(0, 0) * fh.q[0] * p0_u / fh.p[0];
  p_joint(0, 1) = PT(0, 1) * fh.q[1] * p0_u / fh.p[0];
  p_joint(1, 0) = PT(1, 0) * fh.q[0] * (1 - p0_u) / fh.p[1];
  p_joint(1, 1) = PT(1, 1) * fh.q[1] * (1 - p0_u) / fh.p[1];

  const double Z = p_joint(0, 0) + p_joint(0, 1) + p_joint(1, 0) + p_joint(1, 1);
  assert(Z > 0);

  for (size_t u = 0; u < 2; u++)
    for (size_t v = 0; v < 2; v++)
      p_joint(u, v) /= Z;
}


static void
weighted_J_D_branch(const vector<double> &rates, const double T,
                    const FelsHelper &fh,
                    const double p0_parent, double &p0_child,
                    vector<double> &J, vector<double> &D) {

  two_by_two P;
  continuous_time_trans_prob_mat(rates[0], rates[1], T, P);

  // get joint distribution of end states
  two_by_two p_joint;
  joint_post(p0_parent, fh, P, p_joint);
  p0_child = p_joint(0, 0) + p_joint(1, 0);

  // compute weighted expectation of sufficient stats
  two_by_two J0;
  two_by_two J1;
  two_by_two D0;
  two_by_two D1;
  expectation_J(rates[0], rates[1], T, J0, J1);
  expectation_D(rates[0], rates[1], T, D0, D1);

  for (size_t u = 0; u < 2; u++) {
    for (size_t v = 0; v < 2; v++) {
      J[0] += p_joint(u, v)*J0(u, v);
      J[1] += p_joint(u, v)*J1(u, v);
      D[0] += p_joint(u, v)*D0(u, v);
      D[1] += p_joint(u, v)*D1(u, v);
    }
  }
}



/* Iterates over all nodes in pre-order, and for each branch it computes
 * conditional means of sufficient statistics.
 */
static void
downward_process(const TreeHelper &th, const size_t site_id,
                 const vector<vector<Path> > &paths,
                 const vector<double> &init_pi, const vector<double> &rates,
                 const vector<FelsHelper> &fh,
                 vector<vector<double> > &J, vector<vector<double> > &D,
                 vector<double> &init_pi_post) {

  vector<double> p0_margin(th.n_nodes);

  // marginalize at root node
  p0_margin[0] = root_post_prob0(init_pi, fh[0].q);
  init_pi_post[0] += p0_margin[0];
  init_pi_post[1] += 1 - init_pi[0];

  for (size_t node_id = 1; node_id < th.n_nodes; ++node_id)
    weighted_J_D_branch(rates, paths[site_id][node_id].tot_time, fh[node_id],
                        p0_margin[th.parent_ids[node_id]],
                        p0_margin[node_id], J[node_id], D[node_id]);
}


/* Iterates over all nodes in pre-order, and for each branch it samples
 * a new path.
 */
static void
sampling_downward(const TreeHelper &th,const vector<Path> &path,
                  const vector<double> &init_pi, const vector<double> &rates,
                  const vector<FelsHelper> &fh,
                  std::mt19937 &gen, vector<Path> &sampled_path) {

  vector<double> p0_margin(th.n_nodes);
  std::uniform_real_distribution<double> unif(0.0, 1.0);

  // marginalize at root node
  p0_margin[0] = root_post_prob0(init_pi, fh[0].q);
  sampled_path = vector<Path>(th.n_nodes);
  sampled_path.front() = Path(path[1].init_state, 0.0);
  //sampled_path.front() = Path(unif(gen) > p0_margin[0], 0.0);

  for (size_t node_id = 1; node_id < th.n_nodes; ++node_id) {
    const size_t start_state = sampled_path[th.parent_ids[node_id]].end_state();
    const double T = path[node_id].tot_time;

    sampled_path[node_id].init_state = start_state;
    sampled_path[node_id].tot_time = T;

    two_by_two P;
    continuous_time_trans_prob_mat(rates[0], rates[1], T, P);
    const double p0 =
    P[start_state][0] * fh[node_id].q[0] / fh[node_id].p[start_state];
    const size_t sampled_state = (unif(gen) > p0);

    const TwoStateCTMarkovModel ctmm(rates[0], rates[1]);
    assert(end_cond_sample_forward_rejection(ctmm, start_state, sampled_state,
                                             T, gen,
                                             sampled_path[node_id].jumps, 0));
  }
}


////////////////////////////////////////////////////////////////////////////////
///////////////          site-independent E-step       /////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* Obtain conditional means of sufficient statistics */
void
expectation_sufficient_statistics(const vector<double> &rates,
                                  const vector<double> &init_pi,
                                  const TreeHelper &th,
                                  const vector<vector<Path> > &paths,
                                  vector<vector<double> > &J,
                                  vector<vector<double> > &D,
                                  vector<double> &init_pi_post) {

  assert(!paths.empty() && paths[0].size() > 1);
  init_pi_post = {0, 0};

  /* sufficient statistics: J - jump counts, D - spent time */
  J = vector<vector<double> > (th.n_nodes, vector<double> (2, 0.0));
  D = vector<vector<double> > (th.n_nodes, vector<double> (2, 0.0));

  for (size_t site_id = 0; site_id < paths.size(); site_id++) {
    vector<FelsHelper> fh(th.n_nodes);
    upward_process(th, site_id, paths, rates, fh);
    downward_process(th, site_id, paths, init_pi, rates, fh, J, D,
                     init_pi_post);
  }

  init_pi_post[0] /= paths.size();
  init_pi_post[1] /= paths.size();
}


/* Sample new evolution history */
void
sample_paths(const vector<double> &rates, const vector<double> &init_pi,
             const TreeHelper &th,
             const vector<vector<Path> > &paths,
             std::mt19937 &gen, vector<vector<Path> > &sampled_paths) {

  assert(!paths.empty() && paths[0].size() > 1);
  const size_t n_sites = paths.size();
  
  sampled_paths.resize(n_sites);

  for (size_t site_id = 0; site_id < n_sites; site_id++) {
    sampled_paths[site_id].resize(paths[site_id].size());
    vector<FelsHelper> fh(th.n_nodes);
    upward_process(th, site_id, paths, rates, fh);

    vector<Path> site_path;
    sampling_downward(th, paths[site_id], init_pi, rates, fh, gen, site_path);

    for (size_t node_id = 0; node_id < th.n_nodes; ++node_id)
      sampled_paths[site_id][node_id] = site_path[node_id];
  }
}

////////////////////////////////////////////////////////////////////////////////
///////////////          site-independent M-step       /////////////////////////
////////////////////////////////////////////////////////////////////////////////

void
compute_sufficient_statistics(const vector<vector<Path> > &paths,
                              vector<vector<double> > &J,
                              vector<vector<double> > &D) {

  assert(!paths.empty() && paths[0].size() > 1);
  const size_t n_sites = paths.size();
  const size_t n_nodes = paths[0].size();

  J.resize(n_nodes);
  D.resize(n_nodes);

  for (size_t b = 1; b < n_nodes; b++) {

    J[b].resize(2);
    D[b].resize(2);
    vector<double> J_branch = {0.0, 0.0};
    vector<double> D_branch = {0.0, 0.0};

    for (size_t site_id = 0; site_id < n_sites; site_id++) {

      size_t prev_state = paths[site_id][b].init_state;
      double time = 0.0;

      for (size_t j = 0; j < paths[site_id][b].jumps.size(); j++) {
        J_branch[prev_state] += 1;
        D_branch[prev_state] += ( paths[site_id][b].jumps[j] - time );
        prev_state = 1 - prev_state;
        time = paths[site_id][b].jumps[j];
      }
      D_branch[prev_state] += (paths[site_id][b].tot_time - time);
    }

    // get averages
    J[b][0] = J_branch[0] / n_sites;
    J[b][1] = J_branch[1] / n_sites;
    D[b][0] = D_branch[0] / n_sites;
    D[b][1] = D_branch[1] / n_sites;
  }
}


void
estimate_rates(const vector<vector<double> > &J,
               const vector<vector<double> > &D,
               vector<double> &rates, TreeHelper &th) {

  vector<double> J_sum(2, 0.0);
  vector<double> D_sum(2, 0.0);

  for (size_t b = 1; b < th.n_nodes; ++b) {
    J_sum[0] += J[b][0];
    J_sum[1] += J[b][1];

    D_sum[0] += D[b][0];
    D_sum[1] += D[b][1];
  }

  // update rates
  if (D_sum[0] > 0)
    rates[0] = std::max(J_sum[0] / D_sum[0], 10e-6);
  if (D_sum[1] > 0)
    rates[1] = std::max(J_sum[1] / D_sum[1], 10e-6);
}


double
indep_rate_scaling_factor(const vector<double> &rates) {

  // stationary rates pi from T
  vector<double> pi (2, 0.0);
  pi[0] = rates[1] / (rates[0] + rates[1]);
  pi[1] = rates[0] / (rates[0] + rates[1]);

  return pi[0] * rates[0] + pi[1] * rates[1];
}



void
estimate_rates_and_branches(const vector<vector<double> > &J,
                            const vector<vector<double> > &D,
                            vector<double> &rates, TreeHelper &th,
                            vector<vector<Path> > &paths) {

  // estimate rates
  estimate_rates(J, D, rates, th);

  // update branch lengths
  for (size_t b = 1; b < th.n_nodes; ++b)
    th.branches[b] *= (J[b][0] + J[b][1]) / (D[b][0]*rates[0] + D[b][1]*rates[1]);


  // set one change per site per unit time
  const double scale_factor = indep_rate_scaling_factor(rates);
  transform(th.branches.begin(), th.branches.end(), th.branches.begin(),
            std::bind(std::multiplies<double>(), std::placeholders::_1,
                      scale_factor));

  transform(rates.begin(), rates.end(), rates.begin(),
            std::bind(std::divides<double>(), std::placeholders::_1,
                      scale_factor));

  // scale jump times
  for (size_t b = 1; b < paths[0].size(); ++b) {
    for (size_t i = 0; i < paths.size(); ++i) {
      const double scale = th.branches[b] / paths[i][b].tot_time;
      for (size_t j = 0; j < paths[i][b].jumps.size(); ++j)
        paths[i][b].jumps[j] *= scale;
      paths[i][b].tot_time = th.branches[b];
    }
  }
}


void
estimate_root_distribution(const vector<vector<Path> > &paths,
                           vector<double> &init_pi) {

  assert(!paths.empty() && paths[0].size() > 1);

  init_pi.resize(2);

  double N1 = 0.0;
  for(size_t site_id = 0; site_id < paths.size(); site_id++)
    N1 += paths[site_id][1].init_state;

  init_pi[1] = N1 / paths.size();
  init_pi[0] = 1 - init_pi[1];
}
