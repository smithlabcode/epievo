
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
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


static void
get_horizontal_summary_stats(const vector<bool> &seq,
                             vector<vector<double> > &stat) {
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


static string
horiz_summary_str(const vector<bool> &seq) {
  vector<vector<double> > stat;
  get_horizontal_summary_stats(seq, stat);

  std::ostringstream oss;
  oss << '('
      << stat[0][0] << ',' << '\t' << stat[0][1] << ',' << '\t'
      << stat[1][0] << ',' << '\t' << stat[1][1] << ')';
  return oss.str();
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

static void
sum_triplet(const vector<bool> &seq, vector<size_t> &triplet_stat) {
  triplet_stat = vector<size_t>(8, 0);
  for (size_t i = 1; i < seq.size()-1; ++i)
    ++triplet_stat[ord(seq[i-1], seq[i], seq[i+1])];
}

static void
get_expo_rate(const model_param &p, vector<double> &triplet_rate) {
  triplet_rate = vector<double>(8, 0.0);
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      for (size_t k =0; k < 2; ++k)
        triplet_rate[4*i+2*j+k] = exp(p.stationary_logfac[i][1-j] +
                                      p.stationary_logfac[1-j][k] +
                                      p.stationary_logbaseline[i][k]);
}

static void
flip_update(const size_t position, vector<bool> &seq,
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

static void
first_jump(const vector<double> &triplet_rate, std::mt19937 &gen,
           vector<bool> &seq, vector<size_t> &triplet_stat,
           vector<Path> &paths, double &time) {
  const double rate = std::inner_product(triplet_stat.begin(), triplet_stat.end(),
                                         triplet_rate.begin(), 0);
  std::exponential_distribution<double> exp_distr(rate);

  // sample time of first jump
  const double holding_time = exp_distr(gen);

  // sample which state context to jump from
  vector<double> context_prob(triplet_stat.size(), 0.0);
  for (size_t i = 0; i < triplet_stat.size(); ++i)
    context_prob[i] = triplet_stat[i]*triplet_rate[i]/rate;

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
  for (size_t i = 1; !found && i < seq.size() - 1; ++i) {
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
    paths[position].jumps.push_back(time + holding_time);

  // update the current time
  time += holding_time;
}


// use new data structure PatSeq
void
first_jump(const vector<double> &triplet_rate, std::mt19937 &gen,
           PatSeq &patseq, vector<Path> &paths, double &time) {

  const size_t p_size = 8;
  vector<size_t> triplet_stat(p_size, 0);
  for (size_t i = 0; i < p_size; ++i)
    triplet_stat[i] = patseq.get_context_freq(i);

  /* sample time of first jump */
  const double rate = std::inner_product(triplet_stat.begin(), triplet_stat.end(),
                                         triplet_rate.begin(), 0);
  std::exponential_distribution<double> exp_distr(rate);
  const double holding_time = exp_distr(gen);

  if (time + holding_time < paths[0].tot_time) {
    /* sample which state context to jump from */
    vector<double> context_prob(triplet_stat.size(), 0.0);
    for (size_t i = 0; i < p_size; ++i)
      context_prob[i] = triplet_stat[i]*triplet_rate[i]/rate;

    std::discrete_distribution<size_t> multinom(context_prob.begin(),
                                                context_prob.end());
    const size_t context = multinom(gen);
    size_t position = patseq.random_mutate(context, gen);  /*make the jump*/
    paths[position].jumps.push_back(time + holding_time);  /* update paths */
  }

  time += holding_time; /* update time */
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


struct watch_info {
  double node_id;
  double time_val;
  double n_dom;
  double mean_dom_size;
  double dom_size_sd;
  double fraction;
  vector<size_t> patfreq;
};


std::ostream &
operator<<(std::ostream &os, const watch_info &w) {
  os << w.node_id << '\t'
     << w.time_val << '\t'
     << w.n_dom << '\t'
     << w.mean_dom_size << '\t'
     << w.dom_size_sd << '\t'
     << w.fraction << '\t';
  for (size_t ct = 0; ct < 8; ++ct)
    os << w.patfreq[ct] << '\t';
  return os;
}


static void
update_watch_stats(const model_param &p,
                   const PatSeq &patseq,
                   const size_t node_id,
                   const double watch,
                   const double time,
                   vector<watch_info> &w) {

  w.push_back(watch_info());

  vector<size_t> dom_size; // domain sizes
  patseq.to_domain_sizes(dom_size);
  const size_t n_domains = dom_size.size();

  const double dom_tot =
    std::accumulate(dom_size.begin(), dom_size.end(), 0.0);

  const double mds = dom_tot/n_domains; // mean domain size

  // standard deviation
  const double sq_sum =
    std::inner_product(dom_size.begin(), dom_size.end(), dom_size.begin(), 0.0);
  const double stdev = std::sqrt(sq_sum/n_domains - mds*mds);

  w.back().n_dom = n_domains;
  w.back().mean_dom_size = mds;
  w.back().dom_size_sd = stdev;
  w.back().fraction = dom_tot/p.n_site;

  patseq.get_all_context_freq(w.back().patfreq);

  w.back().node_id = node_id;
  w.back().time_val = floor(time/watch)*watch;
}


template <class T> string
mat_tostring(const vector<vector<T> > &m) {
  std::ostringstream oss;
  for (size_t i = 0; i < m.size(); ++i) {
    oss << "[" << std::setw(10) << std::right << m[i][0];
    for (size_t j = 1; j < m[i].size(); ++j)
      oss << ',' << std::setw(10) << std::right << m[i][j];
    oss << "]" << endl;
  }
  return oss.str();
}

template <class T> string
vec_tostring(const vector<T> &v) {
  std::ostringstream oss;
  oss << '[';
  if (!v.empty()) {
    oss << v[0];
    for (size_t i = 1; i < v.size(); ++i)
      oss << ',' << v[i];
  }
  oss << ']';
  return oss.str();
}


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

    model_param p;
    p.read_param(param_file);
    if (VERBOSE)
      cerr << "initial distribution factors:" << endl
           << p.t << endl << mat_tostring(p.init_logfac) << endl;

    vector<vector<double> > init_T;
    convert_parameter(p.init_logfac, init_T);
    if (VERBOSE)
      cerr << "initial sequence Markov chain transition matrix:" << endl
           << mat_tostring(init_T) << endl;

    vector<vector<double> > T;
    convert_parameter(p.stationary_logfac, T);
    if (VERBOSE)
      cerr << "Markov chain transition stationary distribution:" << endl
           << mat_tostring(T) << endl;

    if (!pathfile.empty())
      write_pathfile_header(pathfile);

    /* standard mersenne_twister_engine seeded with rd()*/
    std::random_device rd;
    std::mt19937 gen(rd());

    /* (2) INITIALIZE THE ROOT SEQUENCE */
    // initial sequence
    vector<bool> root_seq;
    sample_init_state(p.n_site, T, gen, root_seq);
    if (VERBOSE)
      cerr << "initial frquencies\n" << horiz_summary_str(root_seq) << endl;

    /* starting context frequencies */
    vector<size_t> triplet_stat;
    sum_triplet(root_seq, triplet_stat);
    if (VERBOSE)
      cerr << "triplet_stat" << endl
           << vec_tostring(triplet_stat) << endl;

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

    /* READ IN THE PHYLOGENETIC TREE */
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
    vector<watch_info> w; // only used if "watching" requested

    /* ITERATE OVER THE NODES IN THE TREE */
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

      if (EXTRA_VERBOSE)
        cerr << "n_jumps=" << n_jumps << "\ttime=" << time << "\t"
             << horiz_summary_str(evolution[node_id]) << endl;

      /* SAMPLE CHANGES ALONG THE CURRENT BRANCH */
      double prev_time = 0;
      while (time < branches[node_id]) {
        if (OPTION == 1) { /* boolean vector */
          first_jump(triplet_rate, gen, evolution[node_id],
                     triplet_stat, paths, time);
          if (EXTRA_VERBOSE && (n_jumps + 1) % 1000 == 0)
            cerr << "n_jumps=" << n_jumps << "\ttime=" << time << "\t"
                 << horiz_summary_str(evolution[node_id]) << endl;
        }
        else { /* pattern pos array*/
          first_jump(triplet_rate, gen, patseq, paths, time);
          if (EXTRA_VERBOSE && (n_jumps + 1) % 1000 == 0 && time < branches[node_id]) {
            cerr << "n_jumps=" << n_jumps << "\ttime=" << time << "\tPatSeq_context_freq:";
            for (size_t ct = 0; ct < 8; ++ct)
              cerr << patseq.get_context_freq(ct) << "\t";
            cerr << endl;
          }

          /* update HMR stats */
          if (keep_watch_info && time - prev_time > watch) {
            update_watch_stats(p, patseq, node_id, watch, time, w);
            prev_time = floor(time/watch)*watch;
          }
        }
        ++n_jumps;
      }
      patseq.to_seq(evolution[node_id]);

      if (VERBOSE)
        cerr << n_jumps << "\t" << branches[node_id] << "\t"
             << horiz_summary_str(evolution[node_id]) << endl;

      if (!pathfile.empty())
        append_to_pathfile(pathfile, paths, node_id);
    }

    if (!outfile.empty())
      write_output(outfile, subtree_sizes, node_names, p, evolution);

    if (keep_watch_info) {
      watchfile = pathfile + ".stats";
      std::ofstream outstat(watchfile.c_str());
      outstat << "branch" << "\t" << "time" << "\t"
              << "n.domain" << "\t" << "mean.domain" << "\t"
              << "sd.domain" << "\t" << "fraction" << "\t"
              << "pattern000" << "\t" << "pattern001" << "\t"
              << "pattern010" << "\t" << "pattern011" << "\t"
              << "pattern100" << "\t" << "pattern101" << "\t"
              << "pattern110" << "\t" << "pattern111" << endl;
      copy(w.begin(), w.end(), std::ostream_iterator<watch_info>(outstat, "\n"));
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
