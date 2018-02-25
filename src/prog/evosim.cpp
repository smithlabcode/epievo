
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>   /* exp, sqrt, pow */
#include <numeric>  /* std::inner_product */

#include <sys/stat.h>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "PhyloTreePreorder.hpp"
#include "path.hpp"  /* related to Path */
#include "param.hpp" /* model_param */
#include "TripletPattern.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;


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

void get_random_sequence(const size_t N, vector<bool> &s, std::mt19937 &gen) {
  s.resize(N, true);
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  for (size_t i = 0; i < N; ++i)
    if (unif(gen) < 0.5)
      s[i] = false;
}

static void
get_horizontal_summary_stats(const vector<bool> &seq,
                             vector<vector<double> > &stat) {
  stat = vector<vector<double> >(2, vector<double>(2, 0.0));

  /* count the number of each type of consecutive pair */
  const size_t n_pairs = seq.size() - 1;
  for (size_t i = 0; i < n_pairs; ++i)
    ++stat[seq[i]][seq[i+1]];  /* implicit conversion */

  /* divide by the total to get estimates of probabilities */
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      stat[i][j] /= n_pairs;
}

static void
sample_init_state(const size_t n_site,  const vector<vector<double> > &T,
                  std::mt19937 &gen, vector<bool> &seq) {
  seq.resize(n_site, true);
  double prob1 = (1.0 - T[0][0])/(2-T[1][1] - T[0][0]);
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  seq[0] = (unif(gen) < prob1);
  for (size_t i = 1; i < n_site ; ++i) {
    double r = unif(gen);
    double p = seq[i-1]? T[1][1] : T[0][0];
    seq[i] = (r <= p)? seq[i-1] : !seq[i-1];
  }
}

static void
gibbs_sample_init_state(const size_t n_site,
                        const vector<vector<double> > &logfac,
                        std::mt19937 &gen, vector<bool> &seq) {
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  for (size_t i = 1; i < n_site - 1; ++i) {
    const size_t s = static_cast<size_t>(seq[i]);
    const double accept_prob =
      exp(logfac[seq[i-1]][1-s] + logfac[1-s][seq[i+1]] -
          logfac[seq[i-1]][s] - logfac[s][seq[i+1]]);
    if (unif(gen) < accept_prob)
      seq[i] = !seq[i];
  }
  vector<vector<double> > stat;
  get_horizontal_summary_stats(seq, stat);
}


void gibbs_sample_path(const model_param &p,  std::mt19937 &gen,
                       vector<Path> &paths) {
  // sampling excludes two endpoints
  std::uniform_int_distribution<int> uni(1,  paths.size() - 2);

  for (size_t n = 1; n < paths.size() - 1; ++n) {
    size_t i = uni(gen);
    Environment env;
    intersect_paths(paths[i-1], paths[i+1], env);

    // Propose a new path
    Path new_path = paths[i];
    new_path.jumps.clear();

    // track current state and current time
    // while generating the new path
    bool s = paths[i].init_state;
    double t = 0.0;

    vector<std::exponential_distribution<double> > exp_distrs(2);
    for (size_t j = 0; j < env.breaks.size(); ++j) {
      for (size_t k = 0; k < 2; ++k) {
        const double rate =
          exp(p.stationary_logfac[env.left[j]][1-k] +
              p.stationary_logfac[1-k][env.right[j]] +
              p.stationary_logbaseline[env.left[j]][env.right[j]]);
        exp_distrs[k] = std::exponential_distribution<double>(rate);
      }
      while (t < env.breaks[j]) {
        const double interval = exp_distrs[s](gen);
        if ((t + interval) <= env.breaks[j]) {
          new_path.jumps.push_back(t + interval);
          s = !s;
          t += interval;
        } else {
          t = env.breaks[j];
        }
      }
    }

    // take care of the last interval
    for (size_t k = 0; k < 2; ++k) {
      assert(env.left.size() && env.right.size());
      const double rate =
        exp(p.stationary_logfac[env.left.back()][1-k] +
            p.stationary_logfac[1-k][env.right.back()]+
            p.stationary_logbaseline[env.left.back()][env.right.back()]);
      exp_distrs[k] = std::exponential_distribution<double>(rate);
    }
    while (t < env.tot_time) {
      const double interval = exp_distrs[s](gen);
      if ((t + interval) < env.tot_time) {
        new_path.jumps.push_back(t + interval);
        s = !s;
        t += interval;
      } else {
        t = env.tot_time;
      }
    }

    // update path
    paths[i] = new_path; // maybe paths[i].swap(new_path);?
  }
}


////////////////////////////////////////////////////////////////////////////////
// Continuous time Markov chain of sequences
////////////////////////////////////////////////////////////////////////////////
size_t ord(bool i, bool j, bool k) { // needs more intuitive function name
  return i*4+j*2+k;
}

void ord_to_state(size_t context,
                  bool &i, bool &j, bool &k) {
  assert(context < 8);
  i = (bool)(context/4);
  context %= 4;
  j = (bool)(context/2);
  context %= 2;
  k = (bool)(context);
}

void sum_triplet(const vector<bool> &seq, vector<size_t> &triplet_stat) {
  triplet_stat = vector<size_t>(8, 0);
  for (size_t i = 1; i < seq.size()-1; ++i)
    ++triplet_stat[ord(seq[i-1], seq[i], seq[i+1])];
}

void get_expo_rate(const model_param &p, vector<double> &triplet_rate) {
  triplet_rate = vector<double>(8, 0.0);
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      for (size_t k =0; k < 2; ++k)
        triplet_rate[4*i+2*j+k] = exp(p.stationary_logfac[i][1-j] +
                                      p.stationary_logfac[1-j][k] +
                                      p.stationary_logbaseline[i][k]);
}

void flip_update(const size_t position, vector<bool> &seq,
                 vector<size_t> &triplet_stat) {
  if (position > 1) {
    --triplet_stat[ord(seq[position-2], seq[position-1], seq[position])];
    ++triplet_stat[ord(seq[position-2], seq[position-1], !seq[position])];
  }
  --triplet_stat[ord(seq[position-1], seq[position], seq[position+1])];
  ++triplet_stat[ord(seq[position-1], !seq[position], seq[position+1])];
  if (position < seq.size()-2) {
    --triplet_stat[ord(seq[position], seq[position+1], seq[position+2])];
    ++triplet_stat[ord(!seq[position], seq[position+1], seq[position+2])];
  }
  seq[position] = !seq[position];
}

void first_jump(const vector<double> &triplet_rate, std::mt19937 &gen,
                vector<bool> &seq, vector<size_t> &triplet_stat,
                vector<Path> &paths, double &time) {
  const double rate = std::inner_product(triplet_stat.begin(), triplet_stat.end(),
                                         triplet_rate.begin(), 0);
  std::exponential_distribution<double> exp_distr(rate);

  /* sample time of first jump */
  const double holding_time = exp_distr(gen);

  /* sample which state context to jump from */
  vector<double> context_prob(triplet_stat.size(), 0.0);
  for (size_t i = 0; i < triplet_stat.size(); ++i) {
    context_prob[i] = triplet_stat[i]*triplet_rate[i]/rate;
  }
  std::discrete_distribution<size_t> multinom(context_prob.begin(),
                                              context_prob.end());
  size_t context = multinom(gen);

  /* sample which position to flip */
  std::uniform_real_distribution<double> unif(0.0,1.0);
  size_t position = 0;
  const size_t n =
    (static_cast<size_t>(unif(gen)*triplet_stat[context]) %
     triplet_stat[context]) + 1;
  bool found = false;
  size_t count = 0;
  bool I, J, K;
  ord_to_state(context, I, J, K);
  for (size_t i = 1; !found && i < seq.size()-1; ++i) {
    // ADS should define some inline functions for this below
    if ( !( (seq[i-1]^I) | (seq[i]^J) | (seq[i+1]^K) ) )
      ++count;

    if (count == n) {
      found = true;
      position = i;
    }
  }
  assert(found);

  // flip the position and update triplet_stat
  flip_update(position, seq, triplet_stat);

  // update paths
  if (time + holding_time < paths[position].tot_time)
    paths[position].jumps.push_back(time+holding_time);
  //update time
  time += holding_time;
}


// use new data structure PatSeq
void first_jump(const vector<double> &triplet_rate,  std::mt19937 &gen,
                PatSeq &patseq, vector<Path> &paths, double &time) {

  const size_t p_size = 8;
  vector<size_t> triplet_stat(p_size, 0);
  for (size_t i = 0; i < p_size; ++i) {
    triplet_stat[i] = patseq.get_context_freq(i);
  }

  /* sample time of first jump */
  const double rate = std::inner_product(triplet_stat.begin(), triplet_stat.end(),
                                         triplet_rate.begin(), 0);
  std::exponential_distribution<double> exp_distr(rate);
  const double holding_time = exp_distr(gen);

  if (time + holding_time < paths[0].tot_time) {
    /* sample which state context to jump from */
    vector<double> context_prob(triplet_stat.size(), 0.0);
    for (size_t i = 0; i < p_size; ++i) {
      context_prob[i] = triplet_stat[i]*triplet_rate[i]/rate;
    }

    std::discrete_distribution<size_t> multinom(context_prob.begin(),
                                                context_prob.end());
    size_t context = multinom(gen);
    size_t position = patseq.random_mutate(context, gen);  /*make the jump*/
    paths[position].jumps.push_back(time + holding_time);  /* update paths */
  }

  time += holding_time; /* update time */
}


////////////////////////////////////////////////////////////////////////////////
// stationary and Markov chain
////////////////////////////////////////////////////////////////////////////////
void convert_parameter(const vector<vector<double> > &stationary_logfac,
                      vector<vector<double> > &T) {
  vector<vector<double> > Q = stationary_logfac;
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      Q[i][j] = exp(stationary_logfac[i][j]);
  const double delta = sqrt(pow(Q[0][0] - Q[1][1], 2) + 4*Q[0][1]*Q[1][0]);
  T = Q; /* transition probability matrix */

  // compute the diagonal entries
  const double diag_denom = Q[0][0]+Q[1][1] + delta;
  T[1][1] = 2*Q[1][1]/diag_denom;
  T[0][0] = 2*Q[0][0]/diag_denom;

  // now compute the anti-diagonal entries
  const double anti_numer = 4*Q[0][1]*Q[1][0];
  T[0][1] = anti_numer/(pow(Q[0][0] + delta, 2) - Q[1][1]*Q[1][1]);
  T[1][0] = anti_numer/(pow(Q[1][1] + delta, 2) - Q[0][0]*Q[0][0]);
}

static
void write_pathfile_header(const string &pathfile) {
  std::ofstream outpath(pathfile.c_str());
  outpath << "## paths" << endl;
}

static void
append_to_pathfile(const string &pathfile, vector<Path> &paths,
                   const size_t node_id) {
  std::ofstream outpath(pathfile.c_str(), std::ofstream::app);
  for (size_t i = 0; i < paths.size(); ++i) {
    outpath << node_id << "\t" <<  i << "\t"
            << paths[i].init_state << "\t" << 0;
    for (size_t j = 0; j < paths[i].jumps.size(); ++j)
      outpath << "," << paths[i].jumps[j];
    if (paths[i].jumps.size() == 0 ||
        paths[i].jumps.back() < paths[i].tot_time)
      outpath << "," << paths[i].tot_time << endl;
    else
      outpath << endl;
  }
}


static void
write_output(const string &outfile, const vector<size_t> &subtree_sizes,
             const vector<string> &node_names, const model_param &p,
             const vector<vector<bool> > &evolution) {
  std::ofstream out(outfile.c_str());
  out << "pos\t";
  for (size_t i = 0; i < subtree_sizes.size(); ++ i)
    out << node_names[i] << "\t";
  out << endl;
  for (size_t pos = 0; pos < p.n_site; ++pos) {
    out << pos ;
    for (size_t i = 0; i < subtree_sizes.size(); ++ i) {
      out << "\t" << (size_t)(evolution[i][pos]);
    }
    out << endl;
  }
}

static void
print_mat(vector<vector<double> > &m) {
  cerr << "[" << endl;
  for (size_t i = 0; i < m.size(); ++i) {
    for (size_t j = 0; j < m[0].size(); ++j)
      cerr << m[i][j] << "\t";
    cerr << endl;
  }
  cerr << "]"<< endl;

}

template <class T> void
print_vec(vector<T> &v) {
  cerr << "[";
  for (size_t i = 0; i < v.size(); ++i)
    cerr << v[i] << ",\t";
  cerr << endl;
}

//////////////////////////////
//////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// SIMULATION
////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char **argv) {

  try {

    string outfile;
    string pathfile;
    bool VERBOSE = false;
    size_t OPTION = 2;
    double watch = 0;
    string watchfile;

    OptionParser opt_parse(strip_path(argv[0]), "simulate methylome evolution",
                           "<params-file>");
    opt_parse.add_opt("output", 'o', "name of output file for methylomes"
                      "(default: stdout)", false, outfile);
    opt_parse.add_opt("paths", 'p', "name of output file for evolution paths"
                      "(default: stdout)", false, pathfile);
    opt_parse.add_opt("watch", 'w', "print summary statistics "
                      "at specified time interval (when -o)", false, watch);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
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

    std::ofstream outstat;
    if (!pathfile.empty() && watch > 0) {
      watchfile = pathfile;
      watchfile.append(".stats");
      outstat.open(watchfile.c_str(), std::ofstream::out);
      outstat << "branch" << "\t" << "time" << "\t"
              << "n.domain" << "\t" << "mean.domain" << "\t"
              << "sd.domain" << "\t" << "fraction" << "\t"
              << "pattern000" << "\t" << "pattern001" << "\t"
              << "pattern010" << "\t" << "pattern011" << "\t"
              << "pattern100" << "\t" << "pattern101" << "\t"
              << "pattern110" << "\t" << "pattern111" << endl;
    }

    /* (1) INITIALIZING PARAMETERS */
    if (VERBOSE)
      cerr << "Reading parameter from file " << param_file << endl;

    model_param p;
    p.read_param(param_file);
    if (VERBOSE)
      cerr << "read in params initial distribution factors" << endl
           << p.t << endl
           << p.init_logfac[0][0] << "\t" << p.init_logfac[0][1] << endl
           << p.init_logfac[1][0] << "\t" << p.init_logfac[1][1] << endl;

    vector<vector<double> > init_T;
    convert_parameter(p.init_logfac, init_T);
    if (VERBOSE) {
      cerr << "Initial sequence Markov chain transition matrix:" << endl;
      print_mat(init_T);
    }

    vector<vector<double> > T;
    convert_parameter(p.stationary_logfac, T);
    if (VERBOSE) {
      cerr << "Stationary sequence Markov chain transition matrix:" << endl;
      print_mat(T);
    }

    if (!pathfile.empty())
      write_pathfile_header(pathfile);

    /*Standard mersenne_twister_engine seeded with rd()*/
    std::random_device rd;
    std::mt19937 gen(rd());

    /* (2) INITIALIZE THE ROOT SEQUENCE */
    /* initial sequence by two-state Markov chain */
    vector<bool> root_seq;
    sample_init_state(p.n_site, T, gen, root_seq);

    /* starting context frequencies */
    vector<size_t> triplet_stat;
    sum_triplet(root_seq, triplet_stat);
    if (VERBOSE) {
      cerr << "triplet_stat:\t";
      print_vec(triplet_stat);
    }

    /* evolution rates (constant throughout the evolution) */
    vector<double> triplet_rate;
    get_expo_rate(p, triplet_rate);

    if (VERBOSE) {
      double tot_rate = 0;
      for (size_t i = 0; i < triplet_rate.size(); ++i) {
        bool I, J, K;
        ord_to_state(i, I, J, K);
        cerr << int(I) << int(J) << int(K) << "\t"
             << triplet_rate[i] << "\t"
             << (double)(triplet_stat[i])/root_seq.size() <<endl;
        tot_rate += triplet_rate[i]*triplet_stat[i]/root_seq.size();
      }
      cerr << "mutation per site (at root): " << tot_rate << endl;
    }

    vector<size_t> subtree_sizes;
    p.t.get_subtree_sizes(subtree_sizes);
    p.t.assign_missing_node_names();
    vector<string> node_names;
    p.t.get_node_names(node_names);
    const size_t n_nodes = subtree_sizes.size();
    vector<size_t> parent_ids;
    get_parent_id(subtree_sizes, parent_ids);
    vector<double> branches;
    p.t.get_branch_lengths(branches);

    if (VERBOSE)
      cerr << "[tree:]\n" << p.t.tostring() << endl;

    vector<vector<bool> > evolution(subtree_sizes.size(), root_seq);
    vector<WatchStat> watch_stats;

    for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
      if (VERBOSE)
        cerr << "node " << node_names[node_id]
             << "\t" << branches[node_id] << endl;

      evolution[node_id] = evolution[parent_ids[node_id]];
      sum_triplet(evolution[node_id], triplet_stat);
      PatSeq patseq(evolution[node_id]);
      double time = 0;
      size_t n_jumps = 0;

      /* initialize paths */
      vector<Path> paths;
      initialize_paths(evolution[node_id], branches[node_id], paths);

      if (VERBOSE) {
        vector<vector<double> > stat;
        get_horizontal_summary_stats(evolution[node_id], stat);
        cerr << n_jumps << "\t" << time << "\t"
             << stat[0][0] << "\t" << stat[0][1] << "\t"
             << stat[1][0] << "\t" << stat[1][1] << endl;
      }

      double prev_time = 0;
      while (time < branches[node_id]) {
        if (OPTION == 1) { /* boolean vector */
          first_jump(triplet_rate, gen, evolution[node_id],
                     triplet_stat, paths, time);
          if (VERBOSE && n_jumps % 1000 == 0) {
            vector<vector<double> > stat;
            get_horizontal_summary_stats(evolution[node_id], stat);
            cerr << n_jumps << "\t" << time << "\t";
            print_mat(stat);
          }
        } else { /* pattern pos array*/
          first_jump(triplet_rate, gen, patseq, paths, time);
          if (VERBOSE && n_jumps%1000 == 0) {
            cerr << n_jumps << "\t" << time << "\t";
            vector<size_t> all_freq;
            patseq.get_all_context_freq(all_freq);
            print_vec(all_freq);
          }

          /* output HMR stat*/
          if (watch > 0 && time - prev_time > watch) {
            double round_time = floor(time/watch)*watch;
            WatchStat ws;
            ws.set(patseq, node_id, round_time, p.n_site);
            watch_stats.push_back(ws);
            prev_time = round_time;
          }
        }
        ++n_jumps;
      }
      patseq.to_seq(evolution[node_id]);

      if (VERBOSE) {
        vector<vector<double> > stat;
        get_horizontal_summary_stats(evolution[node_id], stat);
        cerr << n_jumps << "\t" << branches[node_id] << "\t";
        print_mat(stat);
      }

      if (!pathfile.empty())
        append_to_pathfile(pathfile, paths, node_id);
    }

    if (!outfile.empty())
      write_output(outfile, subtree_sizes, node_names, p, evolution);

    if (!pathfile.empty() && watch > 0) {
      for (size_t i = 0; i <  watch_stats.size(); ++i)
        watch_stats[i].write(outstat);
    }
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
