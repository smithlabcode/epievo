
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <math.h>   /* exp, sqrt, pow */
#include <numeric>  /* std::inner_product */

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "PhyloTreePreorder.hpp"
#include "path.hpp"  /* related to Path */
#include "param.hpp" /* model_param */

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;

bool file_exist(string param_file) {
  std::ifstream in(param_file.c_str());
  return in.good();
}

void get_random_sequence(const size_t N, vector<bool>&s) {
  s.resize(N, true);
  std::random_device rd;
  //Standard mersenne_twister_engine seeded with rd()
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> unif(0.0,1.0);
  for (size_t i = 0; i < N; ++i)
    if (unif(gen) < 0.5)
      s[i] =false;
}

void summary(const vector<bool> &seq, vector<vector<double> > &stat) {
  stat = vector<vector<double> >(2, vector<double>(2, 0.0));
  for (size_t i = 0; i < seq.size()-1; ++i) {
    stat[seq[i]][seq[i+1]] += 1;  //implicit conversion
  }
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      stat[i][j] = stat[i][j]/(seq.size()-1);
}

void gibbs_sample_init_state(const size_t n_site,
                             const vector<vector<double> > &logfac,
                             vector<bool>& seq) {
  std::random_device rd;
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<double> unif(0.0,1.0);
  for (size_t i = 1; i < n_site - 1; ++i) {
    const size_t s = static_cast<size_t>(seq[i]);
    const double accept_prob =
      exp(logfac[seq[i-1]][1-s] + logfac[1-s][seq[i+1]] -
          logfac[seq[i-1]][s] - logfac[s][seq[i+1]]);
    if (unif(gen) < accept_prob)
      seq[i] = !seq[i];
  }
  vector<vector<double> > stat;
  summary(seq, stat);
}


void gibbs_sample_path(const model_param &p, vector<Path>& paths) {
  std::random_device rd;
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<int> uni(1,  paths.size()-2);
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
    double t = 0;

    vector<std::exponential_distribution<double> > exp_distrs(2);
    for (size_t j = 0; j < env.breaks.size(); ++j) {
      for (size_t k = 0; k < 2; ++k) {
        // exponential distribution parameters (add scaler later)
        double rate = exp(p.stationary_logfac[env.left[j]][1-k] +
                          p.stationary_logfac[1-k][env.right[j]] +
                          p.stationary_logbaseline[env.left[j]][env.right[j]]);
        exp_distrs[k] = std::exponential_distribution<double>(rate);
      }
      while (t < env.breaks[j]) {
        double interval = exp_distrs[s](gen);
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
      double rate =
        exp(p.stationary_logfac[env.left.back()][1-k] +
            p.stationary_logfac[1-k][env.right.back()]+
            p.stationary_logbaseline[env.left.back()][env.right.back()]);
      exp_distrs[k] = std::exponential_distribution<double>(rate);
    }
    while (t < env.tot_time) {
      double interval = exp_distrs[s](gen);
      if ((t + interval) < env.tot_time) {
        new_path.jumps.push_back(t + interval);
        s = !s;
        t += interval;
      } else {
        t = env.tot_time;
      }
    }

    // update path
    paths[i] = new_path;
  }
}

////////////////////////////////////////////////////////////////////////////////
// Continuous time Markov chain of sequences
////////////////////////////////////////////////////////////////////////////////
size_t ord(bool i, bool j, bool k) {
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
  std::uniform_real_distribution<double> unif(0.0,1.0);
  const double sample = unif(gen);
  double cdf = context_prob[0];
  size_t context = 0;
  while (sample > cdf && context < 8) {
    ++context;
    cdf += context_prob[context];
  }

  assert (context < 8);

  // sample which position to flip
  size_t position = 0;
  const size_t n = (((size_t)(unif(gen)*triplet_stat[context])) %
                    triplet_stat[context]) + 1;
  bool found = false;
  size_t count = 0;
  bool I, J, K;
  ord_to_state(context, I, J, K);
  for (size_t i = 1; !found && i < seq.size()-1; ++i) {
    if ( !( (seq[i-1]^I) | (seq[i]^J) | (seq[i+1]^K) ) )
      ++count;

    if (count == n) {
      found = true;
      position = i;
    }
  }
  assert(found);

  // flip the position and update triplet_stat
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
  if (time+holding_time < paths[position].tot_time)
    paths[position].jumps.push_back(time+holding_time);
  //update time
  time += holding_time;
}


////////////////////////////////////////////////////////////////////////////////
// stationary and Markov chain
////////////////////////////////////////////////////////////////////////////////
void covert_parameter(const vector<vector<double> > &stationary_logfac,
                      vector<vector<double> > &T) {
  vector<vector<double> > Q = stationary_logfac;
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      Q[i][j] = exp(stationary_logfac[i][j]);
  double delta = sqrt(pow(Q[0][0] - Q[1][1], 2) + 4*Q[0][1]*Q[1][0]);
  // transition probability matrix
  T = Q;
  T[1][1] = 2*Q[1][1] /(Q[0][0]+Q[1][1] + delta );
  T[0][0] = 2*Q[0][0] /(Q[0][0]+Q[1][1] + delta );
  T[0][1] = 4*Q[0][1]*Q[1][0]/(pow(Q[0][0] + delta, 2) -Q[1][1]*Q[1][1]);
  T[1][0] = 4*Q[0][1]*Q[1][0]/(pow(Q[1][1] + delta, 2) -Q[0][0]*Q[0][0]);

}

////////////////////////////////////////////////////////////////////////////////
// SIMULATION
////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char **argv) {

  try {

    string outfile;
    string pathfile;
    bool VERBOSE = false;

    OptionParser opt_parse(strip_path(argv[0]), "simulate methylome evolution",
                           "<params-file>");
    opt_parse.add_opt("output", 'o', "name of output file for methylomes"
                      "(default: stdout)", false, outfile);
    opt_parse.add_opt("paths", 'p', "name of output file for evolution paths"
                      "(default: stdout)", false, pathfile);
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
    if (!file_exist(param_file)) {
      cerr << "File "<< param_file << "doesn't exist." << endl;
      cerr << opt_parse.help_message() << endl;
      return  EXIT_SUCCESS;
    }

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
    covert_parameter(p.stationary_logfac, T);

    if (VERBOSE)
      cerr << "Markov chain transition matrix corresponds "
           << "to stationary distribution" << endl
           << "[" << T[0][0] << "\t" << T[0][1] << endl
           << " " << T[1][0] << "\t" << T[1][1] << "]"<< endl;

    std::ofstream outpath;
    if (!pathfile.empty()){
      outpath.open (pathfile.c_str(), std::ofstream::out);
      outpath << "## paths" << endl;
      outpath.close();
    }

    // initial sequence
    vector<bool> root_seq;
    get_random_sequence(p.n_site, root_seq);
    for (size_t i = 0; i < 500; ++i) {
      gibbs_sample_init_state(p.n_site, p.init_logfac, root_seq);
    }

    if (VERBOSE) {
      vector<vector<double> > stat;
      summary(root_seq, stat);
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
      summary(target_seq, stat);
      cerr << "target frquencies" << endl
           << "(" << stat[0][0] << ",\t" << stat[0][1] << ",\t"
           << stat[1][0] << ",\t"  << stat[1][1] << ")" << endl;
    }

    // starting context frequencies
    vector<size_t> triplet_stat;
    sum_triplet(root_seq, triplet_stat);
    if (VERBOSE) {
      cerr << "triplet_stat" << endl;
      for (size_t i =0; i < triplet_stat.size(); ++i)
        cerr << i << "\t" << triplet_stat[i] << endl;
    }

    // evolution rates (constant throughout the evolution)
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

    for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
      if (VERBOSE)
        cerr << "node " << node_names[node_id]
             << "\t" << branches[node_id] << endl;

      evolution[node_id] = evolution[parent_ids[node_id]];
      sum_triplet(evolution[node_id], triplet_stat);

      double time = 0;
      size_t n_jumps = 0;

      // initialize paths
      vector<Path> paths;
      initialize_paths(evolution[node_id], branches[node_id], paths);

      if (VERBOSE) {
        vector<vector<double> > stat;
        summary(evolution[node_id], stat);
        cout << n_jumps << "\t" << time << "\t"
             << stat[0][0] << "\t" << stat[0][1] << "\t"
             << stat[1][0] << "\t" << stat[1][1] << endl;
      }

      while (time < branches[node_id]) {
        first_jump(triplet_rate, gen, evolution[node_id], triplet_stat, paths, time);
        ++n_jumps;

        if (VERBOSE && n_jumps % 1000 == 0) {
          vector<vector<double> > stat;
          summary(evolution[node_id], stat);
          cout << n_jumps << "\t" << time << "\t"
               << stat[0][0] << "\t" << stat[0][1] << "\t"
               << stat[1][0] << "\t" << stat[1][1] << endl;
        }
      }

      if (VERBOSE) {
        vector<vector<double> > stat;
        summary(evolution[node_id], stat);
        cout << n_jumps << "\t" << branches[node_id] << "\t"
             << stat[0][0] << "\t" << stat[0][1] << "\t"
             << stat[1][0] << "\t" << stat[1][1] << endl;
      }

      if (!pathfile.empty()) {
        outpath.open (pathfile.c_str(), std::ofstream::app);
        for (size_t i = 0; i < paths.size(); ++i) {
          outpath << node_id << "\t" <<  i << "\t"
                  << paths[i].init_state << "\t" << 0;
          for (size_t j = 0; j < paths[i].jumps.size(); ++j)
            outpath << "," << paths[i].jumps[j];
          if (paths[i].jumps.size() == 0 ||
              paths[i].jumps.back() < paths[i].tot_time)
            outpath << "," <<  paths[i].tot_time << endl;
          else
            outpath << endl;
        }
        outpath.close();
      }
    }

    if (!outfile.empty()) {
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
