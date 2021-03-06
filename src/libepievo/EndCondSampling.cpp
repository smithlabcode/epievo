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
#include <algorithm>   // std::lower_bound,
#include <iostream>
#include <cmath>
#include <random>
#include <functional>

#include "EndCondSampling.hpp"
#include "EpiEvoModel.hpp"
#include "Path.hpp"
#include "epievo_utils.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::string;

using std::function;
using std::bind;
using std::placeholders::_1;
using std::begin;
using std::end;

using std::exponential_distribution;
using std::uniform_real_distribution;
using std::ref;

static const double NUMERICAL_TOLERANCE = 1e-10;
static const size_t max_sample_count = 1e10;


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////
////////  #######  #### #######  #######  #####  ########
////////  ##    ##  ##  ##    ## ##      ##   ##    ##
////////  ##    ##  ##  #######  #####   ##         ##
////////  ##    ##  ##  ##   ##  ##      ##   ##    ##
////////  #######  #### ##    ## #######  #####     ##
////////
////////
////////

/* Probability of no jumps: no jumps has nonzero probability only when
 * both end-points have the same state. We have two versions: one
 * assumes the same state at both end points ("a") and the other will
 * take another end point ("b") and return probabilty 0.0 if a !=
 * b. The expression for Pr(no jumps) is given by:
 *
 * p_a = exp(-Qa*T)/Paa(T)
 *
 * which can be found, for example, in Hobolth & Thorne (2014),
 * equation 1.13.
 */
double
prob_no_jump(const CTMarkovModel &the_model,
             const two_by_two &P, const double T, const size_t a) {
  return exp(-the_model.get_rate(a)*T)/P(a, a);
}

double
prob_no_jump(const CTMarkovModel &the_model,
             const two_by_two &P, const double T,
             const size_t a, const size_t b) {
  return (a == b) ? prob_no_jump(the_model, P, T, a) : 0.0;
}


///////////
///////////  THE NEXT 3 FUNCTIONS EXIST ONLY FOR SANITY CHECK
///////////

/* This function is J_{aj} in Hobolth & Stone (2009), which appears on
 * page 1209, below eqn (2.6). In our case, with only two states, we
 * don't need to handle the case of (lambda_j + Q_a = 0), but an
 * assertion checks anyway.
 */
static double
helper_integral_total_Jaj(const double T,      // T (HS2009)
                          const double lambda, // lambda_j (HS2009)
                          const double Qa) {   // Q_a (HS2009)

  // this assertion checks a condition that should always hold for our
  // case of two states, but for the general case the value of J_{ai}
  // can be evaluated using Te^{T\lambda_j} if this does not hold. But
  // we don't care...
  assert(lambda + Qa != 0.0);

  // there is a good chance some refactoring of the expression below
  // will allow evaluation over a larger domain without numerical
  // problems.
  const double numer = (exp(T*lambda) - exp(-Qa*T));
  assert(std::isfinite(numer));
  const double r = numer/(lambda + Qa);
  assert(std::isfinite(r));

  return r;
}

/* This function computes the value in the summation of eqn (2.6) of
 * Hobolth & Stone (2009). The coefficient of Q_{ai}/P_{ab}(T) is not
 * included here, and the value "j" from that paper may only take
 * values of 0 or 1 for our binary state space.
 */
double
summation_in_total_cdf(const CTMarkovModel &the_model,
                       const double T, const size_t a, const size_t b) {

  const size_t a_bar = complement_state(a);
  const double a_rate = the_model.get_rate(a);

  const double integrand_0 =
    helper_integral_total_Jaj(T, the_model.eigen_values[0], a_rate);
  const double eigen_mat_factor_0 = the_model.U[a_bar][0]*the_model.Uinv[0][b];

  const double integrand_1 =
    helper_integral_total_Jaj(T, the_model.eigen_values[1], a_rate);
  const double eigen_mat_factor_1 = the_model.U[a_bar][1]*the_model.Uinv[1][b];

  const double r = (eigen_mat_factor_0*integrand_0 +
                    eigen_mat_factor_1*integrand_1);

  assert(std::isfinite(r));
  assert(r > 0.0);

  return r;
}

/* this total cumulative density is eqn (2.4) in HS2009, and should
 * integrate to 1.0 if a != b. Otherwise I'm not exactly sure what it
 * is doing, except that it is the largest value the integral from 0
 * to T can take, so any sampled CDF for the distribution must be no
 * larger than this value.
 */
static double
total_cumulative_density(const CTMarkovModel &the_model,
                         const two_by_two &P,
                         const double T, const size_t a, const size_t b) {

  const double Qai = the_model.get_rate(a); // two-states so Qa = Qai
  return (Qai/P(a, b))*summation_in_total_cdf(the_model, T, a, b);
}

///////////
///////////  FUNCTIONS BELOW DO THE WORK FOR DIRECT SAMPLING
///////////

/* This function appears in Hobolth & Stone (2009), page 1209, as part
 * of REMARK 4. In our case, with only two states, we don't need to
 * handle the case of (lambda_j + Q_a = 0). The function is used for
 * the inverse transformation method, where a cdf probability is
 * sampled, and the inverse of the CDF is used to determine a sample
 * value that follows the corresponding PDF.
 *
 * Arguments:
 *
 *  T =      total time for the interval
 *  lambda = eigenvalue corresponding to state "j" in eqn (2.5) of HS2009
 *  Qa =     rate for holding in start state "a" (Qa = -Qaa = Sum Qab for a!=b)
 *  t =      some time less than T for which we want to evaluate the CDF
 *
 * Returns: the value of the integral; not sure what conditions it
 * must satisfy except that it should be maximized when t=T
 */
static double
helper_integral(const double T,
                const double lambda,  // lambda_j (HS2009)
                const double Qa,      // Q_a (HS2009)
                const double t) {     // t (HS2009)

  // condition below is required becuase we divide by this value.
  assert(lambda + Qa != 0.0);

  // some helpers that are used more than once.
  const double neg_lambda_plus_Qa = -1.0*(lambda + Qa);
  const double T_x_lambda = T*lambda;

  /* This "numerator" expression has the sign reversed for the two
   * terms, which requires a correction below. This way seems to allow
   * the expression to remain finite over a larger domain of its
   * arguments.
   */
  const double numer = exp(T_x_lambda + t*neg_lambda_plus_Qa) - exp(T_x_lambda);
  assert(std::isfinite(numer));

  /* Below the denominator has the negative sign, which does not
   * appear in HS2009, but is needed because our numerator above has
   * an extra negative sign.
   */
  const double r = numer/neg_lambda_plus_Qa;
  assert(std::isfinite(r));

  return r;
}


/* This function computes the value in the summation of eqn (2.5) of
 * Hobolth & Stone (2009). The coefficient of Q_{ai}/P_{ab}(T) is not
 * included here, and the value "j" from that paper may only take
 * values of 0 or 1 for our binary state space.
 */
double
summation_in_cdf(const CTMarkovModel &the_model,
                 const double T, const size_t a, const size_t b, const double t) {

  const size_t a_bar = complement_state(a);
  const double Qa = the_model.get_rate(a);

  const double integr_0 = helper_integral(T, the_model.eigen_values[0], Qa, t);
  const double eigen_mat_factor_0 = the_model.U[a_bar][0]*the_model.Uinv[0][b];

  const double integr_1 = helper_integral(T, the_model.eigen_values[1], Qa, t);
  const double eigen_mat_factor_1 = the_model.U[a_bar][1]*the_model.Uinv[1][b];

  const double r = (eigen_mat_factor_0*integr_0 + eigen_mat_factor_1*integr_1);

  assert(std::isfinite(r));
  assert(r > 0.0);

  return r;
}

double
cumulative_density_function(const CTMarkovModel &the_model,
                            const two_by_two &P,
                            const double T, const size_t a, const size_t b,
                            const double t) {

  const double Qai = the_model.get_rate(a); // two-states so Qa = Qai
  return (Qai/P(a, b))*summation_in_cdf(the_model, T, a, b, t);
}

/* The summation in the pdf is given by equation (2.5) in HS2009. The
 * pdf is the function f_i(t) and this function evaluates the
 * summation, which lacks the coefficient Qai/Pab(T). This function is
 * used for evaluating probabilities associated with existing paths.
 */
double
summation_in_pdf(const CTMarkovModel &the_model,
                 const double T, const size_t a, const size_t b,
                 const double t) {

  const size_t a_bar = complement_state(a);
  const double Qa = the_model.get_rate(a);

  const double lambda_0 = the_model.eigen_values[0];
  const double lambda_1 = the_model.eigen_values[1];

  const double factor_0 = the_model.U[a_bar][0]*the_model.Uinv[0][b];
  const double factor_1 = the_model.U[a_bar][1]*the_model.Uinv[1][b];

  const double f =
    (factor_0 * exp((T - t) * lambda_0 - t * Qa)) +
    (factor_1 * exp((T - t) * lambda_1 - t * Qa));

  assert(std::isfinite(f));
  assert(f > 0.0);

  return f;
}


/* PDF, probability density function for the first jump happening at a
 * particular time
 */
double
probability_density_function(const CTMarkovModel &the_model,
                             const two_by_two &P,
                             const double T, const size_t a, const size_t b,
                             const double t) {

  //const double Qai = the_model.get_rate(a); // two-states so Qa = Qai
  const double pi = summation_in_total_cdf(the_model, T, a, b);
  return summation_in_pdf(the_model, T, a, b, t) / pi;
}

static double
bisection_search_cumulative_density(const CTMarkovModel &the_model,
                                    const two_by_two &P,
                                    const double T,
                                    const size_t a, const size_t b,
                                    const double target) {

  double lo = 0.0, lo_val = 0.0; // equals modified cdf for val 0.0

  double hi = T;
  // use the "total" version below just because it's faster
  double hi_val = summation_in_total_cdf(the_model, T, a, b);

  // ADS: the target is a value u sampled from (0, x), where x is the
  // cumulative density for the first change being to state i!=a over
  // the total time T. We use the "summation_in_cdf" here to avoid
  // having to repeatedly multiply the Qai/Pab each time we evaluate
  // any density. This also seems to help keep values centered. At the
  // same time we need to modify the "target" accordingly with the
  // coefficient Qai/Pab.
  const double transformed_target = target/(the_model.get_rate(a)/P[a][b]);

  assert(lo_val <= transformed_target && hi_val >= transformed_target);

  while (hi - lo > NUMERICAL_TOLERANCE) {
    const double mi = (lo + hi)/2.0;
    const double mi_val = summation_in_cdf(the_model, T, a, b, mi);
    if (mi_val >= transformed_target) {
      hi = mi;
      hi_val = mi_val;
    }
    else {
      lo = mi;
      lo_val = mi_val;
    }
  }
  return (lo + hi)/2.0;
}


/* Continuous time Markov chian with rate matrix Q.
   Return the first jump time within (0, T) or T if no jumps,
   given state at time 0 being a, and state at T being b.
*/
double
end_cond_sample_first_jump(const CTMarkovModel &the_model,
                           const size_t a, const size_t b,
                           const double T, function<double()> &unif) {

  assert(T > std::numeric_limits<double>::min());

  two_by_two P; // P = exp(QT)
  the_model.get_trans_prob_mat(T, P);

  // if (a == b) then we decide whether or not to have any jumps
  if (a == b && unif() < prob_no_jump(the_model, P, T, a))
    return T;

  // if we have jumps, we must sample the waiting time tau in state a
  // using inverse transformation method. x~pdf <=> CDF(x)~Unif(0, 1)
  const double cdf_range = total_cumulative_density(the_model, P, T, a, b);

  // should use a root finder here because the bisection seems
  // extremely slow
  const double first_jump_time =
    bisection_search_cumulative_density(the_model, P, T, a, b, unif()*cdf_range);

  return first_jump_time;
}


/* Endpoint-conditioned sampling of path witin time interval T */
void
end_cond_sample_direct(const TwoStateCTMarkovModel &the_model,
                       const size_t start_state, const size_t end_state,
                       const double T, std::mt19937 &gen,
                       vector<double> &jump_times, const double start_time) {

  std::uniform_real_distribution<double> unif(0.0, 1.0);
  function<double()> distr(bind(unif, std::ref(gen)));

  const CTMarkovModel ctmm(the_model.rate0, the_model.rate1);

  size_t current_state = start_state;
  double consumed_time =
    end_cond_sample_first_jump(ctmm, current_state, end_state, T, distr);

  // ADS: the use of NUMERICAL_TOLERANCE below should be checked. We
  // need to make sure that the sampling of a jump will return exactly
  // the total time interval when it should, and not some
  // approximation to it.
  while (T - consumed_time > NUMERICAL_TOLERANCE) {
    jump_times.push_back(start_time + consumed_time);
    current_state = complement_state(current_state);
    consumed_time +=
      end_cond_sample_first_jump(ctmm, current_state,
                                 end_state, T - consumed_time, distr);
  }
}


/* Direct sampling and return proposal densities */
void
end_cond_sample_direct_prop(const TwoStateCTMarkovModel &the_model,
                            const size_t start_state, const size_t end_state,
                            const double T, std::mt19937 &gen,
                            vector<double> &jump_times, double &prob,
                            const double start_time) {

  prob = 1.0;
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  function<double()> distr(bind(unif, std::ref(gen)));

  const CTMarkovModel ctmm(the_model.rate0, the_model.rate1);
  two_by_two PT;
  ctmm.get_trans_prob_mat(T, PT);

  size_t current_state = start_state;
  double rest_time = T;
  double tau =
  end_cond_sample_first_jump(ctmm, current_state, end_state, T, distr);
  double p_no_jump = prob_no_jump(ctmm, PT, rest_time, current_state);
  double consumed_time = tau;

  while (T - consumed_time > NUMERICAL_TOLERANCE) {
    if (current_state == start_state) {
      p_no_jump = prob_no_jump(ctmm, PT, rest_time, current_state);
      prob *= (1 - p_no_jump);
    }
    prob *= probability_density_function(ctmm, PT, rest_time,
                                         current_state, end_state, tau);
    jump_times.push_back(start_time + consumed_time);
    current_state = complement_state(current_state);
    rest_time -= tau;

    tau = end_cond_sample_first_jump(ctmm, current_state, end_state,
                                     rest_time, distr);
    consumed_time += tau;
    if (rest_time > NUMERICAL_TOLERANCE)
      ctmm.get_trans_prob_mat(rest_time, PT);
  }
  prob *= prob_no_jump(ctmm, PT, rest_time, current_state);

}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////
/////// #######  ######  #######  ##     ##    ##    #######  #######
/////// ##      ##    ## ##    ## ##  #  ##   ####   ##    ## ##    ##
/////// #####   ##    ## #######  ##  #  ## ##    ## #######  ##    ##
/////// ##      ##    ## ##   ##  ##  #  ## ######## ##   ##  ##    ##
/////// ##       ######  ##    ##  ### ###  ##    ## ##    ## #######
////////
//////// Hobolth & Stone (2009)

inline double rexp_inv(const double mu, std::mt19937 &gen,
                       uniform_real_distribution<double> &unif) {
  return (- 1 / mu * log(unif(gen)));
}

size_t
forward_sampling(vector<function<double()> > &the_distrs,
                 size_t a, const double T, const double start_time,
                 vector<double> &jump_times) {
  double tau = start_time;
  while ((tau += the_distrs[a]()) < T) {
    a = complement_state(a);
    jump_times.push_back(tau);
  }
  return a;
}


bool
end_cond_sample_forward_rejection(exponential_distribution<double> &exp0,
                                  exponential_distribution<double> &exp1,
                                  const size_t start_state,
                                  const size_t end_state,
                                  const double T,
                                  std::mt19937 &gen,
                                  vector<double> &jump_times,
                                  const double start_time) {
  
  vector<function<double()> > the_distrs = {
    function<double()>(bind(exp0, ref(gen))),
    function<double()>(bind(exp1, ref(gen)))
  };

  size_t sample_count = 0;
  vector<double> proposal;
  while (forward_sampling(the_distrs, start_state, T, 0.0,
                          proposal) != end_state &&
         sample_count < max_sample_count) {
    ++sample_count;
    proposal.clear();
  }

  if (sample_count < max_sample_count)
    transform(begin(proposal), end(proposal), std::back_inserter(jump_times),
              bind(std::plus<double>(), _1, start_time));
  assert((proposal.size() % 2) == static_cast<size_t>(start_state != end_state));

  return (sample_count < max_sample_count);
}


bool
end_cond_sample_forward_rejection(const TwoStateCTMarkovModel &the_model,
                                  const size_t start_state,
                                  const size_t end_state,
                                  const double T,
                                  std::mt19937 &gen,
                                  vector<double> &jump_times,
                                  const double start_time) {

  typedef exponential_distribution<double> exp_distr;
  vector<function<double()> > the_distrs = {
    function<double()>(bind(exp_distr(the_model.rate0), ref(gen))),
    function<double()>(bind(exp_distr(the_model.rate1), ref(gen)))
  };

  size_t sample_count = 0;
  vector<double> proposal;
  while (forward_sampling(the_distrs, start_state, T, 0.0,
                          proposal) != end_state &&
         sample_count < max_sample_count) {
    ++sample_count;
    proposal.clear();
  }

  if (sample_count < max_sample_count)
    transform(begin(proposal), end(proposal), std::back_inserter(jump_times),
              bind(std::plus<double>(), _1, start_time));
  assert((proposal.size() % 2) == static_cast<size_t>(start_state != end_state));

  return (sample_count < max_sample_count);
}


bool
end_cond_sample_forward_rejection(const TwoStateCTMarkovModel &the_model,
                                  const size_t start_state,
                                  const size_t end_state,
                                  const double T,
                                  std::mt19937 &gen,
                                  uniform_real_distribution<double> &unif,
                                  vector<double> &jump_times,
                                  const double start_time) {

  vector<function<double()> > the_distrs = {
    function<double()>(bind(rexp_inv, the_model.rate0, ref(gen),ref(unif))),
    function<double()>(bind(rexp_inv, the_model.rate1, ref(gen), ref(unif)))
  };
  size_t sample_count = 0;
  vector<double> proposal;
  while (forward_sampling(the_distrs, start_state, T, 0.0,
                          proposal) != end_state &&
         sample_count < max_sample_count) {
    ++sample_count;
    proposal.clear();
  }
  
  if (sample_count < max_sample_count)
    transform(begin(proposal), end(proposal), std::back_inserter(jump_times),
              bind(std::plus<double>(), _1, start_time));
  assert((proposal.size() % 2) == static_cast<size_t>(start_state != end_state));
  
  return (sample_count < max_sample_count);
}

/* Used for inverse transform sampling, from Nielsen (2001), eqn (A2) */
static double
sample_trunc_exp(function<double()> &U, const double lambda, const double T) {
  return -log(1.0 - U()*(1.0 - exp(-lambda*T)))/lambda;
}


bool
end_cond_sampling_Nielsen(const TwoStateCTMarkovModel &the_model,
                          size_t start_state, const size_t end_state,
                          const double T, std::mt19937 &gen,
                          vector<double> &jump_times, const double start_time) {

  // if start and end state are the same, use regular forward sampling
  if (start_state == end_state)
    return end_cond_sample_forward_rejection(the_model, start_state, end_state,
                                             T, gen, jump_times, start_time);

  // otherwise do the modified procedure
  typedef std::exponential_distribution<double> exp_distr;
  vector<function<double()> > distr = {
    function<double()>(bind(exp_distr(the_model.rate0), ref(gen))),
    function<double()>(bind(exp_distr(the_model.rate1), ref(gen)))
  };
  function<double()> U(bind(uniform_real_distribution<double>(0, 1), ref(gen)));

  size_t sample_count = 0;
  vector<double> proposal;
  bool valid_path = false;
  while (!valid_path && sample_count < max_sample_count) {
    proposal.clear();
    size_t prev_state = start_state;
    const double Qa = the_model.get_rate(prev_state);
    const double offset = sample_trunc_exp(U, Qa, T);
    proposal.push_back(offset);
    prev_state = complement_state(prev_state);
    if (forward_sampling(distr, prev_state, T, offset, proposal) == end_state)
      valid_path = true;
    ++sample_count;
  }

  if (proposal.size() > 0)
    transform(begin(proposal), end(proposal), std::back_inserter(jump_times),
              bind(std::plus<double>(), _1, start_time));

  return valid_path;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// # #  ###  ###  ###   #   ##   # #  ###  ###   #   ###  ###   #   ###
// # #  # #   #   #    # #  # #  ###   #     #  # #   #    #   # #  # #
// # #  # #   #   ##   # #  ##   ###   #    #   ###   #    #   # #  # #
// # #  # #   #   #    # #  # #  # #   #   #    # #   #    #   # #  # #
// ###  # #  ###  #     #   # #  # #  ###  ###  # #   #   ###   #   # #

struct CTMarkovUnif {
  // CTMarkovUnif is not aware of epigenomic histories
  CTMarkovUnif(const TwoStateCTMarkovModel &the_model);
  double unif_trans_prob(const size_t state_a, const size_t state_b,
                         const size_t n) const;
  size_t us;
  double scaler;
  double r;
};


CTMarkovUnif::CTMarkovUnif(const TwoStateCTMarkovModel &the_model) {
  us = (the_model.rate0 < the_model.rate1);
  scaler = the_model.get_rate(us);
  r = the_model.get_rate(!us) / the_model.get_rate(us);
}

double
CTMarkovUnif::unif_trans_prob(const size_t state_a, const size_t state_b,
                              const size_t n) const {
  const double r_sign = (n % 2 == 0) ? 1.0 : -1.0;
  double prob_stay;
  if (state_a == us)
    prob_stay = (r + r_sign * pow(r, n)) / (1 + r);
  else
    prob_stay = (1 + r_sign * pow(r, n+1)) / (1 + r);
  return (state_a == state_b) ? prob_stay : 1.0 - prob_stay;
}

static size_t
num_unif_trans(const TwoStateCTMarkovModel &the_model,
               const CTMarkovUnif &unif_model,
               const double T, const size_t state_a, const size_t state_b,
               function<double()> &unif, double &prob) {
  two_by_two P; // P = exp(QT)
  continuous_time_trans_prob_mat(the_model.rate0, the_model.rate1, T, P);

  const double u = unif();

  const double muT = unif_model.scaler * T;
  const double nom_const = (state_b == unif_model.us) ? unif_model.r : 1.0;
  const double nom_sign = (state_b == unif_model.us) ? 1.0 : -1.0;
  double nom_series = (state_a == unif_model.us) ? 1 : - unif_model.r;
  const double denom = 1 + unif_model.r;

  size_t n = 0;
  double prob_pois = exp(- muT) / P(state_a, state_b);
  double prob_unif = 1;

  prob = prob_pois * prob_unif * (state_a == state_b);
  double sum_probs = prob;

  while (sum_probs < u) {
    ++n;
    prob_pois *= ( muT / n);
    nom_series *= - unif_model.r;
    prob_unif = (nom_const + nom_sign * nom_series) / denom;
    prob = prob_pois * prob_unif;
    sum_probs += prob;
  }
  return n;
}


void
end_cond_sample_unif(const TwoStateCTMarkovModel &the_model,
                     const size_t start_state, const size_t end_state,
                     const double T, std::mt19937 &gen,
                     vector<double> &jump_times, const double start_time) {

  std::uniform_real_distribution<double> unif(0.0, 1.0);
  function<double()> distr(bind(unif, std::ref(gen)));

  CTMarkovUnif unif_model(the_model);

  // sample number of transitions
  double prob_n_trans;

  size_t num_trans = num_unif_trans(the_model, unif_model, T,
                                    start_state, end_state, distr, prob_n_trans);

  size_t prev_state = start_state;

  double prob = prob_n_trans;
  if (num_trans < 1) {
    // no jumps
    assert(start_state == end_state);
  }
  else if (num_trans == 1) {
    // one jump
    prob *= 1 / T;
    const double trans_time = start_time + distr() * T;
    if (start_state != end_state) {
      // real jump
      jump_times.push_back(trans_time);
      prev_state = complement_state(prev_state);
    }
  }
  else {
    // at least two jumps
    vector<double> trans_times; // all jumps, including virtual jumps
    for (size_t i = 0; i < num_trans; i++) {
      trans_times.push_back(distr() * T);
    }
    std::sort(trans_times.begin(), trans_times.end());

    // determine the state of jumps
    size_t n_real_jumps = 0;
    for (size_t i = 0; i < num_trans - 1; i++) {
      // PDF of arriving time part
      prob *= (i + 1) / T;

      const double next_end = unif_model.unif_trans_prob(!prev_state, end_state,
                                                         num_trans - i - 1);
      const double prev_end = unif_model.unif_trans_prob(prev_state, end_state,
                                                         num_trans - i);
      const double prob_jump = unif_model.unif_trans_prob(prev_state,
                                                          !prev_state, 1) *
      next_end / prev_end;
      if (distr() < prob_jump) {
        // true jump sampled
        prev_state = complement_state(prev_state);
        jump_times.push_back(start_time + trans_times[i]);
        prob *= prob_jump;
        n_real_jumps++;
      } else {
        prob *= (1 - prob_jump);
      }
    }
    // last jump
    if (prev_state != end_state) {
      prev_state = complement_state(prev_state);
      jump_times.push_back(start_time + trans_times.back());
      n_real_jumps++;
    }
  }
  assert(prev_state == end_state);
}


void
end_cond_sample_unif_prop(const TwoStateCTMarkovModel &the_model,
                          const size_t start_state, const size_t end_state,
                          const double T, std::mt19937 &gen,
                          vector<double> &jump_times, double &prob,
                          const double start_time) {

  std::uniform_real_distribution<double> unif(0.0, 1.0);
  function<double()> distr(bind(unif, std::ref(gen)));

  CTMarkovUnif unif_model(the_model);

  // sample number of transitions
  double prob_n_trans;

  size_t num_trans = num_unif_trans(the_model, unif_model, T,
                                    start_state, end_state, distr, prob_n_trans);

  size_t prev_state = start_state;

  prob = prob_n_trans;
  if (num_trans < 1) {
    // no jumps
    assert(start_state == end_state);
  } else if (num_trans == 1) {
    // one jump
    prob *= 1 / T;
    const double trans_time = start_time + distr() * T;
    if (start_state != end_state) {
      // real jump
      jump_times.push_back(trans_time);
      prev_state = complement_state(prev_state);
    }
  } else {
    // at least two jumps
    vector<double> trans_times; // all jumps, including virtual jumps
    for (size_t i = 0; i < num_trans; i++) {
      trans_times.push_back(distr() * T);
    }
    std::sort(trans_times.begin(), trans_times.end());

    // determine the state of jumps
    size_t n_real_jumps = 0;
    for (size_t i = 0; i < num_trans - 1; i++) {
      // PDF of arriving time part
      prob *= (i + 1) / T;

      const double next_end = unif_model.unif_trans_prob(!prev_state, end_state,
                                                         num_trans - i - 1);
      const double prev_end = unif_model.unif_trans_prob(prev_state, end_state,
                                                         num_trans - i);
      const double prob_jump = unif_model.unif_trans_prob(prev_state,
                                                          !prev_state, 1) *
      next_end / prev_end;
      if (distr() < prob_jump) {
        // true jump sampled
        prev_state = complement_state(prev_state);
        jump_times.push_back(start_time + trans_times[i]);
        prob *= prob_jump;
        n_real_jumps++;
      } else {
        prob *= (1 - prob_jump);
      }
    }
    // last jump
    if (prev_state != end_state) {
      prev_state = complement_state(prev_state);
      jump_times.push_back(start_time + trans_times.back());
      n_real_jumps++;
    }
  }
  assert(prev_state == end_state);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////
/////// ########   #######  ####  ######   ######   #######  ##    ##
/////// ##     ## ##     ##  ##  ##       ##       ##     ## ####  ##
/////// ########  ##     ##  ##   ######   ######  ##     ## ## ## ##
/////// ##        ##     ##  ##        ##       ## ##     ## ##  ####
/////// ##         #######  ####  ######   ######   #######  ##    ##

double
expected_num_jumps(const TwoStateCTMarkovModel &the_model,
                   const size_t start_state, const size_t end_state,
                   const double T) {

  const double r0 = the_model.rate0;
  const double r1 = the_model.rate1;
  const double s = r0 + r1;
  const double p = r0 * r1;
  const double d = r1 - r0;
  const double e = exp(- s * T);
  double N = 0;

  if (start_state == end_state) {
    N = 2 * p / s;
    if (start_state == 0)
      N *= ( ( (r1 - r0 * e) * T  - d * (1 - e) / s) / (r1 + r0 * e) );
    else
      N *= ( ( (r0 - r1 * e) * T  + d * (1 - e) / s) / (r0 + r1 * e) );
  } else
    N = 2 * p * T * (1 + e) / (s * (1 - e)) + d * d / (s * s);

  return (N > 0 ? N : (s * T / 2));
}


static size_t
num_Poisson_trans(const double rate, const double T, const size_t state_a,
                  const size_t state_b, function<double()> &unif, double &prob) {
  const double u = unif();

  const double muT = rate * T;
  const double denom = (state_a == state_b) ?
  (exp(muT) + exp(-muT)) / 2 : (exp(muT) - exp(-muT)) / 2;

  size_t n = (state_a == state_b) ? 0 : 1;
  prob = (state_a == state_b) ? 1 : muT;
  double sum_probs = prob;

  while (sum_probs < u * denom) {
    n += 2;
    prob *= ( muT * muT / (n * (n-1)));
    sum_probs += prob;
  }
  return n;
}


void
end_cond_sample_Poisson(const TwoStateCTMarkovModel &the_model,
                        const size_t start_state, const size_t end_state,
                        const double T, std::mt19937 &gen,
                        vector<double> &jump_times, const double start_time) {

  const double rate = expected_num_jumps(the_model, start_state, end_state, T)/T;
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  function<double()> distr(bind(unif, std::ref(gen)));

  // sample number of transitions
  double prob_n_trans;
  size_t num_trans = num_Poisson_trans(rate, T, start_state, end_state, distr,
                                       prob_n_trans);
  assert(((start_state == end_state) && (num_trans % 2 == 0)) ||
         ((start_state != end_state) && (num_trans % 2 != 0)));

  if (num_trans > 0) {
    // sample the time points of jumps
    vector<double> trans_times;
    for (size_t i = 0; i < num_trans; i++)
      trans_times.push_back(distr() * T);
    std::sort(trans_times.begin(), trans_times.end());

    // update the jumps vector
    for (size_t i = 0; i < num_trans; i++)
      jump_times.push_back(start_time + trans_times[i]);
  }
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// #####  #####   ####  #####   ####   ###   #   ##     #####  #####   ####  #####
// ##  ## ##  ## ##  ## ##  ## ##  ## #     ###  ##     ##  ## ##  ## ##  ## ##  ##
// #####  #####  ##  ## #####  ##  ## #### #   # ##     #####  #####  ##  ## #####
// ##     ## ##  ##  ## ##     ##  ##    # ##### ##     ##     ## ##  ##  ## ##  ##
// ##     ##  ##  ####  ##      ####  ###  #   # #####  ##     ##  ##  ####  #####

/* (log) probability density. */
double
end_cond_sample_prob(const TwoStateCTMarkovModel &the_model,
                     const vector<double> &jump_times,
                     const size_t start_state, const size_t end_state,
                     const double start_time, const double end_time,
                     const size_t start_jump, const size_t end_jump) {

  // start_jump must specify the first jump after the start time
  assert(start_jump == 0 || jump_times[start_jump - 1] < start_time);
  assert(start_jump == jump_times.size() || jump_times[start_jump] > start_time);

  // end_jump must specify the first jump after the end time
  assert(end_jump == 0 || jump_times[end_jump - 1] < end_time);
  assert(end_jump == jump_times.size() || jump_times[end_jump] > end_time);

  // the end points should be equal of we have an even number of jumps
  // inside the interval
  assert((end_jump - start_jump) % 2 == static_cast<size_t>(start_state != end_state));

  double log_p = 0.0;

  double curr_time = start_time;

  // if start_jump == end_jump then no jump exists within the
  // specified time interval; otherwise start_jump must specify a time
  // inside the interval
  size_t a = start_state;
  const vector<double> rates = {the_model.rate0, the_model.rate1};
  const vector<double> log_rates = {log(the_model.rate0), log(the_model.rate1)};

  for (size_t i = start_jump; i < end_jump; ++i) {

    const double tau = jump_times[i] - curr_time;
    const double jump_log_prob = log_rates[a] - rates[a] * tau;
    log_p += jump_log_prob;

    assert(std::isfinite(log_p));
    a = complement_state(a);
    curr_time = jump_times[i];
  }

  // calculate last interval
  const double tau = end_time - curr_time;
  assert(a == end_state);
  // last interval
  log_p += - rates[a] * tau;
  log_p -= log(the_model.get_trans_prob(end_time - start_time,
                                        start_state, end_state));
  return log_p;
}
