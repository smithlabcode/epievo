
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <math.h>       /* exp */

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::string;

struct model_param {
  size_t n_site;
  double tot_time;
  vector<vector<double> > stationary_logfac;
  vector<vector<double> > stationary_logbaseline;
  vector<vector<double> > init_logfac;
};


void read_param(const string param_file, model_param &p) {
  std::ifstream in(param_file.c_str());
  string dummy_label;
  in >> dummy_label >> p.n_site;
  in >> dummy_label >> p.tot_time;
  p.stationary_logfac =
    vector<vector<double> >(2, vector<double>(2, 0.0));
  p.stationary_logbaseline =
    vector<vector<double> >(2, vector<double>(2, 0.0));
  p.init_logfac =
    vector<vector<double> >(2, vector<double>(2, 0.0));
  in >> dummy_label;
  in >> p.stationary_logfac[0][0] >> p.stationary_logfac[0][1];
  in >> p.stationary_logfac[1][0] >> p.stationary_logfac[1][1];
  in >> dummy_label;
  in >> p.stationary_logbaseline[0][0] >> p.stationary_logbaseline[0][1];
  in >> p.stationary_logbaseline[1][0] >> p.stationary_logbaseline[1][1];
  in >> dummy_label;
  in >> p.init_logfac[0][0] >> p.init_logfac[0][1];
  in >> p.init_logfac[1][0] >> p.init_logfac[1][1];
}


void get_random_sequence(const size_t N, vector<bool>&s) {
  s.resize(N, true);
  std::random_device rd;
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
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
  // cerr << "(" << stat[0][0] << ",\t" << stat[0][1] << ",\t"
  //      << stat[1][0] << ",\t"  << stat[1][1] << ")" << endl;
}

////////////////////////////////////////////////////////////////////////////////
// path
////////////////////////////////////////////////////////////////////////////////

struct path{
  bool init_state;
  double tot_time;
  vector<double> jumps;
};

struct environment{
  vector<bool> left;
  vector<bool> right;
  vector<double> breaks;
  double tot_time;
};
  
void initialize_paths(const model_param &p, const vector<bool> &seq,
                      vector<path> &paths) {
  paths.resize(seq.size());
  for (size_t i = 0; i < seq.size(); ++i) {
    paths[i].init_state = seq[i];
    paths[i].tot_time = p.tot_time;
    paths[i].jumps.resize(0);
  }
}

void intersect_paths(const path &pa, const path &pb, environment &env) {
  assert (pa.tot_time == pb.tot_time);
  bool sa = pa.init_state;
  bool sb = pb.init_state;
  size_t i = 0;
  size_t j = 0;
  env.tot_time = pa.tot_time;
  while (i < pa.jumps.size() || j < pb.jumps.size()) {  
    env.left.push_back(sa);
    env.right.push_back(sb);
    if (i < pa.jumps.size() && j < pb.jumps.size()) {
      if (pa.jumps[i] < pb.jumps[j]) {
        env.breaks.push_back(pa.jumps[i]);
        ++i;
        sa = !sa; 
      } else if (pa.jumps[i] > pb.jumps[j]) {
        env.breaks.push_back(pb.jumps[j]);
        ++j;
        sb = !sb; 
      } else {
        env.breaks.push_back(pb.jumps[j]);
        ++j;
        ++i;
        sa = !sa;
        sb = !sb;
      }
    } else if (i < pa.jumps.size()) {
      env.breaks.push_back(pa.jumps[i]);
      ++i;
      sa = !sa;
    } else {
      env.breaks.push_back(pb.jumps[j]);
      ++j;
      sb = !sb;
    }
  }
  
  if (env.breaks.size() == 0 ||
      env.breaks.back() < env.tot_time) {
    env.left.push_back(sa); 
    env.right.push_back(sb); 
  }
}


void gibbs_sample_path(const model_param &p, vector<path>& paths) {
  std::random_device rd;
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  
  for (size_t i = 1; i < paths.size() - 1; ++i) {
    environment env;    
    intersect_paths(paths[i-1], paths[i+1], env);

    // Propose a new path
    path new_path = paths[i];
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
    // if (paths[i].jumps.size()) {
    //   for (size_t n = 0; n < paths[i].jumps.size(); ++n)
    //     cerr << paths[i].jumps[n] << ";";
    //   cerr << endl;
    // }
  }
}


void end_sequence(const vector<path> &paths,
                  vector<bool> &seq) {
  seq.resize(paths.size());
  for (size_t i = 0; i < paths.size(); ++i) {
    bool s = paths[i].init_state;
    seq[i] = (paths[i].jumps.size() % 2 == 0)? s : !s;     
  }
}


bool state_at_time(const path &p, const double t) {
  bool s = p.init_state;
  for (size_t i = 0; i < p.jumps.size() && p.jumps[i] < t; ++i)
    s = !s;
  return s;
}

void sequence_at_time(const vector<path> &paths,
                      const double t,
                      vector<bool> &seq) {
  seq.resize(paths.size());
  for (size_t i = 0; i < paths.size(); ++i) {
    seq[i] = state_at_time(paths[i], t);
  }
}

int main(int argc, const char **argv){

  try{

    string outfile;
    bool VERBOSE = false;
  
    OptionParser opt_parse(strip_path(argv[0]), "simulate methylome evolution",
                            "<params-file>");
    opt_parse.add_opt("output", 'o', "name of output file "
                      "(default: stdout)", false, outfile);
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

    model_param p;
    read_param(param_file, p);
    cerr << "read in params initial distribution factors" << endl;
    cerr << p.init_logfac[0][0] << "\t" << p.init_logfac[0][1] << endl
         << p.init_logfac[1][0] << "\t" << p.init_logfac[1][1] << endl;

    // initial sequence
    vector<bool> seq;
    get_random_sequence(p.n_site, seq); 
    for (size_t i = 0; i < 500; ++i) {
      gibbs_sample_init_state(p.n_site, p.init_logfac, seq);
    }
    vector<vector<double> > stat;
    summary(seq, stat);
    cerr << "initial frquencies" << endl
         << "(" << stat[0][0] << ",\t" << stat[0][1] << ",\t"
         << stat[1][0] << ",\t"  << stat[1][1] << ")" << endl;


    //target sequence
    vector<bool> target_seq;
    get_random_sequence(p.n_site, target_seq); 
    for (size_t i = 0; i < 500; ++i) {
      gibbs_sample_init_state(p.n_site, p.stationary_logfac, target_seq);
    }
    summary(target_seq, stat);
    cerr << "target frquencies" << endl
         << "(" << stat[0][0] << ",\t" << stat[0][1] << ",\t"
         << stat[1][0] << ",\t"  << stat[1][1] << ")" << endl;


    // path pa, pb;
    // pa.init_state = true;
    // pb.init_state = false;                         
    // pa.tot_time = 10;
    // pb.tot_time = 10;
    // pa.jumps = vector<double>(3, 0.0);
    // pb.jumps = vector<double>(4, 0.0);
    // pa.jumps[0] = 3.0; pa.jumps[1] = 5.0; pa.jumps[2] = 6.0;
    // pb.jumps[0] = 2.0; pb.jumps[1] = 5.0; pb.jumps[2] = 6.5; pb.jumps[3] = 7.0;
    // environment env; 
    // intersect_paths(pa, pb, env);

    vector<bool> end_seq;
    vector<path> paths;
    initialize_paths(p, seq, paths);
    for (size_t i = 0; i < 1000; ++i) {
      gibbs_sample_path(p, paths);
      end_sequence(paths, end_seq);
      summary(end_seq, stat);
      cerr << "(" << stat[0][0] << ",\t" << stat[0][1] << ",\t"
           << stat[1][0] << ",\t"  << stat[1][1] << ")" << endl;
    }

    if (!outfile.empty()) {
      std::ofstream out(outfile.c_str());
      for (size_t i = 0; i < paths.size(); ++i) {
        out << i << "\t" << paths[i].init_state << "\t" << 0;
        for (size_t j = 0; j < paths[i].jumps.size(); ++j)
          out << "," << paths[i].jumps[j];
        if (paths[i].jumps.size() == 0 ||
            paths[i].jumps.back() < paths[i].tot_time)
          out << "," <<  paths[i].tot_time << endl;
        else
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
