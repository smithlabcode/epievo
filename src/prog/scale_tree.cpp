/* Copyright (C) 2019 University of Southern California
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
#include <exception>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "PhyloTreePreorder.hpp"
#include "TreeHelper.hpp"
#include "EpiEvoModel.hpp"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::runtime_error;


int main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;

    string outfile;
    string param_file;
    double scale_factor = 1.0;

    OptionParser opt_parse(strip_path(argv[0]), "scale phylogenetic tree",
                           "<newick-format>");
    opt_parse.add_opt("factor", 'f', "factor to scale by (default: 1.0)",
                      false, scale_factor);
    opt_parse.add_opt("param", 'p', "input file of epievo model parameters",
                      false, param_file);
    opt_parse.add_opt("output", 'o', "output file (default: stdout)",
                      false, outfile);
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string tree_file(leftover_args.front());
    ////////////////////////////////////////////////////////////////////////

    if (VERBOSE)
      cerr << "[READING TREE: " << tree_file << "]" << endl;
    PhyloTreePreorder the_tree; // tree topology and branch lengths
    std::ifstream tree_in(tree_file.c_str());
    if (!tree_in || !(tree_in >> the_tree))
      throw runtime_error("cannot read tree file: " + tree_file);
    TreeHelper th(the_tree);

    if (!param_file.empty()) {
      if (VERBOSE)
        cerr << "[READING PARAMETER FILE: " << param_file << endl;
      EpiEvoModel the_model;
      read_model(param_file, the_model);
      if (VERBOSE)
        cerr << the_model << endl;
      
      const double rates_scale_factor =
        rate_scaling_factor(the_model.triplet_rates);
      // scale rates
      transform(the_model.triplet_rates.begin(), the_model.triplet_rates.end(),
                the_model.triplet_rates.begin(),
                std::bind(std::divides<double>(), std::placeholders::_1,
                          rates_scale_factor));
      // scale tree
      transform(th.branches.begin(), th.branches.end(), th.branches.begin(),
                std::bind(std::multiplies<double>(), std::placeholders::_1,
                          rates_scale_factor));
    }

    transform(th.branches.begin(), th.branches.end(), th.branches.begin(),
              std::bind(std::multiplies<double>(), std::placeholders::_1,
                        scale_factor));
    
    the_tree.set_branch_lengths(th.branches);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    out << the_tree << endl;
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
}
