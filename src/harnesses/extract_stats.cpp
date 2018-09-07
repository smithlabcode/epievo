/* Copyright (C) 2018 University of Southern California
 *                    Jianghan Qu, Andrew D Smith and Xiaojing Ji
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
#include <bitset>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "Path.hpp"
#include "ParamEstimation.hpp"
#include "StateSeq.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;


int main(int argc, const char **argv) {

  try {

    string outfile;
    bool VERBOSE = false;

    ////////////////////////////////////////////////////////////////////////
    OptionParser opt_parse(strip_path(argv[0]), "extract summary statistics"
                           "from complete data (site-specific paths)",
                           "<path-file>");
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("output", 'o', "output stat file",
                      false, outfile);
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string path_file(leftover_args.front());
    ////////////////////////////////////////////////////////////////////////

    if (VERBOSE)
      cerr << "[READING PATHS FILE: " << path_file << "]" << endl;
    vector<vector<Path> > all_paths; // along multiple branches
    vector<string> node_names;
    read_paths(path_file, node_names, all_paths);

    if (VERBOSE)
      cerr << "[EXTRACTING SUMMARY STATISTICS]" << endl;
      
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    if (!out)
      throw std::runtime_error("bad output file: " + outfile);
    
    vector<double> J;
    vector<double> D;
    get_sufficient_statistics(all_paths, J, D);
    for (size_t i = 0; i < 8; i++) {
      const string triplet = std::bitset<3>(i).to_string();
      out << triplet << "\t" << J[i] << "\t" << D[i] << endl;
    }

  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
