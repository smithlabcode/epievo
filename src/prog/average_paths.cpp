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
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <bitset>
#include <random>
#include <functional>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "Path.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::runtime_error;
using std::numeric_limits;
using std::to_string;

typedef vector<vector<double> > matrix;

static void
add_paths(const vector<vector<Path> > &paths, vector<matrix> &average_paths,
          const size_t n_points) {
  for (size_t b = 1; b < paths.size(); b++) {
    const double bin = paths[b].front().tot_time/(n_points - 1);
    for (size_t site_id = 0; site_id < paths[b].size(); site_id++) {
      size_t prev_state = paths[b][site_id].init_state;
      average_paths[b][site_id].front() += prev_state;
      double curr_time = bin;
      for (size_t i = 1; i < n_points; i++) {
        average_paths[b][site_id][i] +=
        paths[1][site_id].state_at_time(curr_time);
        curr_time += bin;
      }
    }
  }
}


static void
write_output(const string &outfile, const vector<matrix> &states,
             const vector<string> &node_names,
             const vector<double> &branch_len) {
  std::ofstream out(outfile);
  if (!out)
    throw runtime_error("bad output file: " + outfile);
  out << "NODE:" << node_names[0] << endl;
  for (size_t b = 1; b < node_names.size(); b++) {
    out << "NODE:" << node_names[b] << "\t" << branch_len[b] << endl;
     for (size_t site_id = 0; site_id < states[b].size(); ++site_id) {
      out << states[b][site_id].front();
      for (size_t t = 1; t < states[b][site_id].size(); t++)
        out << "\t" << states[b][site_id][t];
      out << endl;
    }
  }
}


int main(int argc, const char **argv) {
  try {
    bool VERBOSE = false;

    string outfile;
    size_t n_points = 100;

    ///////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]), "average local paths",
                           "<input directory>");
    opt_parse.add_opt("outfile", 'o', "output file", true, outfile);
    opt_parse.add_opt("npoints", 'n', "number of bins", false, n_points);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.set_show_defaults();
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
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    const string paths_dir(leftover_args[0]);
    ///////////////////////////////////////////////////////////////////////////

    if (VERBOSE)
      cerr << "[READING PATHS FROM: " << paths_dir << "]" << endl;
    vector<string> files;
    read_dir(paths_dir, "local_paths", files);
    const size_t n_files = files.size();

    vector<vector<Path> > paths; // along multiple branches
    vector<string> node_names;
    vector<matrix> average_paths;
    read_paths(files[0], node_names, paths);

    const size_t n_nodes = paths.size();
    const size_t n_sites = paths[1].size();
    vector<double> branch_len(n_nodes);
    average_paths.resize(n_nodes);
    for (size_t b = 1; b < paths.size(); b++) {
      average_paths[b].resize(n_sites);
      branch_len[b] = paths[b].front().tot_time;
      for (size_t site_id = 0; site_id < paths[b].size(); site_id++)
        average_paths[b][site_id].resize(n_points, 0.0);
    }

    add_paths(paths, average_paths, n_points);

    for (size_t i = 1; i < n_files; i++) {
      paths.clear(); // along multiple branches
      vector<string> node_names;
      read_paths(files[i], node_names, paths);
      add_paths(paths, average_paths, n_points);
    }

    /* AVERAGING PATHS */
    for (size_t b = 1; b < n_nodes; b++)
      for (size_t site_id = 0; site_id < n_sites; ++site_id)
        transform(begin(average_paths[b][site_id]),
                  end(average_paths[b][site_id]),
                  begin(average_paths[b][site_id]),
                  [&](const double x) {return x/n_files;});

    if (VERBOSE)
      cerr << "[WRITING OUTPUT TO: " << outfile << "]" << endl;
    write_output(outfile, average_paths, node_names, branch_len);
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
