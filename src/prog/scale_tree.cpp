/* Copyright (C) 2018 University of Southern California
 *                    Jianghan Qu and Andrew D Smith
 *
 * Author: Andrew D. Smith and Jianghan Qu
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

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "PhyloTreePreorder.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;

int main(int argc, const char **argv) {

  try {

    string outfile;
    double scale_factor = 0.0;

    OptionParser opt_parse(strip_path(argv[0]), "scale phylogenetic tree",
                           "<newick-format>");
    opt_parse.add_opt("factor", 'f', "factor to scale by (default: unit branch sum)",
                      false, scale_factor);
    opt_parse.add_opt("output", 'o', "output file (default: stdout)",
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
    const string tree_file(leftover_args.front());
    ////////////////////////////////////////////////////////////////////////

    PhyloTreePreorder t;

    std::ifstream in(tree_file.c_str());
    if (!in || !(in >> t))
      throw std::runtime_error("bad tree file: " + tree_file);

    vector<double> branches;
    t.get_branch_lengths(branches);

    const double total = (scale_factor > 0.0) ? scale_factor :
      std::accumulate(branches.begin(), branches.end(), 0.0);
    std::transform(branches.begin(), branches.end(), branches.begin(),
                   std::bind2nd(std::divides<double>(), total));

    t.set_branch_lengths(branches);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    out << t << endl;
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
