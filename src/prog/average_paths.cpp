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


static void
add_paths(const vector<vector<Path> > &paths,
          vector<vector<double> > &average_paths) {
  const size_t n_points = average_paths[0].size();
  const double bin = paths[1][0].tot_time/(n_points - 1);
  for (size_t site_id = 0; site_id < paths[1].size(); site_id++) {
    size_t prev_state = paths[1][site_id].init_state;
    average_paths[site_id][0] += prev_state;
    double curr_time = bin;
    for (size_t i = 1; i < average_paths[site_id].size(); i++) {
        average_paths[site_id][i] +=
      paths[1][site_id].state_at_time(curr_time);
        curr_time += bin;
    }
  }
}


static void
write_output(const string &outfile, const vector<vector<double> >  &states) {
  std::ofstream out(outfile);
  if (!out)
    throw runtime_error("bad output file: " + outfile);

  const size_t n_sites = states.size();
  for (size_t site_id = 0; site_id < n_sites; ++site_id) {
    out << states[site_id].front();
    for (size_t t = 1; t < states[site_id].size(); t++)
      out << "\t" << states[site_id][t];
    out << endl;
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
    read_paths(files[0], node_names, paths);

    const size_t n_sites = paths[1].size();
    vector<vector<double> > average_paths(n_sites, vector<double>(n_points, 0));
    add_paths(paths, average_paths);

    for (size_t i = 1; i < n_files; i++) {
      paths.clear(); // along multiple branches
      vector<string> node_names;
      read_paths(files[i], node_names, paths);
      add_paths(paths, average_paths);
    }

    /* AVERAGING PATHS */
    for (size_t site_id = 0; site_id < average_paths.size(); ++site_id)
      transform(begin(average_paths[site_id]), end(average_paths[site_id]),
                begin(average_paths[site_id]),
                [&](const double x) {return x/n_files;});

    if (VERBOSE)
      cerr << "[WRITING OUTPUT TO: ]" << outfile << "]" << endl;
    write_output(outfile, average_paths);
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
