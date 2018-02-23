
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


void get_random_sequence(const size_t N, vector<bool> &s) {
  s.resize(N, true);
  std::random_device rd;
  //Standard mersenne_twister_engine seeded with rd()
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  for (size_t i = 0; i < N; ++i)
    if (unif(gen) < 0.5)
      s[i] = false;
}

static void
get_horizontal_summary_stats(const vector<bool> &seq, vector<vector<double> > &stat) {
  stat = vector<vector<double> >(2, vector<double>(2, 0.0));

  // count the number of each type of consecutive pair
  const size_t n_pairs = seq.size() - 1;
  for (size_t i = 0; i < n_pairs; ++i)
    ++stat[seq[i]][seq[i+1]];  // implicit conversion

  // divide by the total to get estimates of probabilities
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      stat[i][j] /= n_pairs;
}

static void
gibbs_sample_init_state(const size_t n_site,
                        const vector<vector<double> > &logfac,
                        vector<bool> &seq) {

  std::random_device rd;
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
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


void gibbs_sample_path(const model_param &p, vector<Path> &paths) {
  std::random_device rd;
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<int> uni(1,  paths.size() - 2); // why subtract 2?
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
        // exponential distribution parameters (add scaler later)
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

void first_jump(const vector<double> &triplet_rate,  std::mt19937 &gen,
                vector<bool> &seq, vector<size_t> &triplet_stat,
                vector<Path> &paths, double &time) {
  const double rate = std::inner_product(triplet_stat.begin(), triplet_stat.end(),
                                         triplet_rate.begin(), 0);
  std::exponential_distribution<double> exp_distr(rate);

  // sample time of first jump
  const double holding_time = exp_distr(gen);

  // sample which state context to jump from
  vector<double> context_prob(triplet_stat.size(), 0.0);
  for (size_t i = 0; i < triplet_stat.size(); ++i) {
    context_prob[i] = triplet_stat[i]*triplet_rate[i]/rate;
  }

  std::discrete_distribution<size_t> multinom(context_prob.begin(),
                                              context_prob.end());
  size_t context = multinom(gen);

  std::uniform_real_distribution<double> unif(0.0,1.0);
  // const double sample = unif(gen);
  // double cdf = context_prob[0];
  // size_t context = 0;
  // while (sample > cdf && context < 8) {
  //   ++context;
  //   cdf += context_prob[context];
  // }
  // assert (context < 8);

  /* sample which position to flip */
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
  // ADS: seems like the stuff below should be in a separate function
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

    // std::uniform_real_distribution<double> unif(0.0,1.0);
    // const double sample = unif(gen);
    // double cdf = context_prob[0];
    // size_t context = 0;
    // while (sample > cdf && context < p_size) {
    //   ++context;
    //   cdf += context_prob[context];
    // }
    // assert (context < p_size);

    size_t position = patseq.random_mutate(context);  /*make the jump*/
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
  // transition probability matrix
  T = Q;

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
    read_param(param_file, p);

    if (VERBOSE)
      cerr << "read in params initial distribution factors" << endl
           << p.t << endl
           << p.init_logfac[0][0] << "\t" << p.init_logfac[0][1] << endl
           << p.init_logfac[1][0] << "\t" << p.init_logfac[1][1] << endl;

    vector<vector<double> > T;
    convert_parameter(p.stationary_logfac, T);

    if (VERBOSE)
      cerr << "Markov chain transition matrix corresponds "
           << "to stationary distribution" << endl
           << "[" << T[0][0] << "\t" << T[0][1] << endl
           << " " << T[1][0] << "\t" << T[1][1] << "]"<< endl;

    if (!pathfile.empty())
      write_pathfile_header(pathfile);

    /* (2) INITIALIZE THE ROOT SEQUENCE */
    // initial sequence
    vector<bool> root_seq;
    get_random_sequence(p.n_site, root_seq);
    for (size_t i = 0; i < 500; ++i) {
      gibbs_sample_init_state(p.n_site, p.init_logfac, root_seq);
    }

    if (VERBOSE) {
      vector<vector<double> > stat;
      get_horizontal_summary_stats(root_seq, stat);
      cerr << "initial frquencies" << endl
           << "(" << stat[0][0] << ",\t" << stat[0][1] << ",\t"
           << stat[1][0] << ",\t"  << stat[1][1] << ")" << endl;
    }

    //target sequence
    vector<bool> target_seq;
    get_random_sequence(p.n_site, target_seq);
    for (size_t i = 0; i < 500; ++i) {
      gibbs_sample_init_state(p.n_site, p.stationary_logfac, target_seq);
    }

    if (VERBOSE) {
      vector<vector<double> > stat;
      get_horizontal_summary_stats(target_seq, stat);
      cerr << "target frquencies" << endl
           << "(" << stat[0][0] << ",\t" << stat[0][1] << ",\t"
           << stat[1][0] << ",\t"  << stat[1][1] << ")" << endl;
    }

    /* starting context frequencies */
    vector<size_t> triplet_stat;
    sum_triplet(root_seq, triplet_stat);
    if (VERBOSE) {
      cerr << "triplet_stat" << endl;
      for (size_t i =0; i < triplet_stat.size(); ++i)
        cerr << i << "\t" << triplet_stat[i] << endl;
    }

    /* evolution rates (constant throughout the evolution) */
    vector<double> triplet_rate;
    get_expo_rate(p, triplet_rate);

    if (VERBOSE) {
      double tot_rate = 0;
      for (size_t i = 0; i < triplet_rate.size(); ++i) {
        bool I, J, K;
        ord_to_state(i, I, J, K);
        cerr << int(I)<< int(J) << int(K) << "\t"
             << triplet_rate[i] << "\t"
             << (double)(triplet_stat[i])/root_seq.size() <<endl;
        tot_rate += triplet_rate[i]*triplet_stat[i]/root_seq.size();
      }
      cerr << tot_rate << " mutation per site (at root)" << endl;
    }

    std::random_device rd;
    std::mt19937 gen(rd());

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

    vector<double> watch_node;
    vector<double> watch_time;
    vector<double> watch_nds; /*number of domains*/
    vector<double> watch_mds;
    vector<double> watch_ds_stdev;
    vector<double> watch_fraction;
    vector<vector<size_t> > watch_patfreq;

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
            cerr << n_jumps << "\t" << time << "\t"
                 << stat[0][0] << "\t" << stat[0][1] << "\t"
                 << stat[1][0] << "\t" << stat[1][1] << endl;
          }
        } else { /* pattern pos array*/
          first_jump(triplet_rate, gen, patseq, paths, time);
          if (VERBOSE && n_jumps%1000 == 0) {
            cerr << n_jumps << "\t" << time << "\t";
            for (size_t ct = 0; ct < 8; ++ct)
              cerr << patseq.get_context_freq(ct) << "\t";
            cerr << endl;
          }

          /* output HMR stat*/
          if (watch > 0 && time - prev_time > watch) {
            vector<size_t> ds; /*domain sizes*/
            patseq.to_domain_sizes(ds);
            double sum = std::accumulate(ds.begin(), ds.end(), 0.0);
            double mds = sum/ds.size();
            double sq_sum = std::inner_product(ds.begin(), ds.end(), ds.begin(), 0.0);
            double stdev = std::sqrt(sq_sum/ds.size() - mds*mds);

            watch_node.push_back(node_id);
            watch_time.push_back(floor(time/watch)*watch);
            watch_nds.push_back(ds.size());
            watch_mds.push_back(mds);
            watch_ds_stdev.push_back(stdev);
            watch_fraction.push_back(sum/p.n_site );
            vector<size_t> f;
            patseq.get_all_context_freq(f);
            watch_patfreq.push_back(f);
            prev_time = floor(time/watch)*watch;
          }
        }
        ++n_jumps;
      }
      patseq.to_seq(evolution[node_id]);

      if (VERBOSE) {
        vector<vector<double> > stat;
        get_horizontal_summary_stats(evolution[node_id], stat);
        cerr << n_jumps << "\t" << branches[node_id] << "\t"
             << stat[0][0] << "\t" << stat[0][1] << "\t"
             << stat[1][0] << "\t" << stat[1][1] << endl;
      }

      if (!pathfile.empty())
        append_to_pathfile(pathfile, paths, node_id);
    }

    if (!outfile.empty())
      write_output(outfile, subtree_sizes, node_names, p, evolution);

    if (!pathfile.empty() && watch > 0) {
      for (size_t i = 0; i < watch_node.size(); ++i) {
        outstat << watch_node[i] << "\t" << watch_time[i] << "\t"
                << watch_nds[i] << "\t" << watch_mds[i] << "\t"
                << watch_ds_stdev[i] << "\t" << watch_fraction[i] << "\t";
        for (size_t ct = 0; ct < 8; ++ct)
          outstat << watch_patfreq[i][ct] << "\t";
        outstat << endl;
      }
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
