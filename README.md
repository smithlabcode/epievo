Epievo includes tools for simulation and inference of epigenome evolution. 

Building and Installing 
=======================

First, download `epievo` in recursive mode,
```
git clone --recursive https://github.com/andrewdavidsmith/epievo.git
```
then install `epievo` by running
```
cd epievo
make install
```

### Other dependencies
To run Python and R scripts included,
`Python2.x` and `R` need to be installed and added to environment variables.
You will also need to have R package `optparse`,
`ggplot2`
and Python libraries `argparse` installed.

File formats
========================
### Model parameters
Our model parameters are organized in below format:
```
stationary  0.85    0.9
baseline        -0.5    -1.5
init    0.85        0.89
```
The stationary line includes stationary distribution of horizontal Markov transition probabilities T00 and T11. Baseline parameters control the symmetric-context mutation rates r0_0 and r1_1. The third line is similar to the first line, which
refers to horizontal transitions in root epigenome.

### Tree format
`EpiEvo` supports [Newick](http://evolution.genetics.washington.edu/phylip/newicktree.html)
format for tree representation.

### Epigenomic states
`EpiEvo` use below format to present epigenomic state at each (aligned) position in each species.
Species labels are consistent to node labels in the tree file, and get sorted in preorder: 
```
#NODE1    NODE2     NODE3     ...
site1   state_site1_node1   state_site1_node2   state_site1_node3
site2   state_site2_node1   state_site2_node2   state_site2_node3
...
```
The header line has 1 fewer fields than rest lines.
After the header line, the first column shows genomic positions.
The rest columns show binary states in corresponding species and genomic site.

### Global jumps
`EpiEvo` use global jumps to retrieve positions and times of mutations.
Global jumps are sorted by time then position:
```
ROOT:NODE1
[Sequence of binary states]
NODE:NODE2
mut1_time   mut1_site
mut2_time   mut2_site
...
NODE:NODE3
...
```
The block of root node will only include a sequence of binary states.
Then, each node block contains a list of time and position of mutation events.

### Local path
Local path is another way to record mutation events.
Different from global jumps, the local path is a list of mutation times at each position.
The format is below:
```
NODE:NODE1
NODE:NODE2
site1   site1_initial_state   site1_total_time   site1_mut1_time    site1_mut2_time   ...
site2   site2_initial_state   site2_total_time   site2_mut1_time    site2_mut2_time   ...
...
NODE:NODE3
...
```
Again, the root node block has no mutation information.


Usage
========================
### Simulate epigenome evolution
Below command will simulate epigenome evolution from given evolution parameters,
tree topology, or the total time duration of a single branch. Users can also
fix the root sequence by specifying `-r` option. 
```
./bin/epievo_sim [OPTIONS] <params-file>
```
Below is the list of available options:
```
  -o, -output    name of output file for states of each node species (default: stdout) 
  -n, -n-sites   length of sequence to simulate (default: 100) 
  -p, -paths     name of output file for evolution paths as sorted jump times 
                 (default: stdout) 
  -s, -seed      rng seed 
  -r, -root      root states file 
  -t, -tree      Newick format tree file 
  -T, -evo-time  evolutionary time 
  -l, -leaf      write only leaf states 
  -S, -scale     scale the branch lengths 
  -R, -rates     use triplet transition rates instead of model parameters
  -v, -verbose   print more run info 
```

### Estimate model parameters from local paths
Below command will run a gradient-ascent method to obtain maximum-likelihood estimates
of model parameters and branch lengths, from provided local paths of complete evolution history:
```
./bin/epievo_sim [OPTIONS] <path-file>
```
Below is the list of available options:
```
  -p, -param            initial parameter file 
  -t, -tree             initial tree file in newick format 
  -v, -verbose          print more run info 
  -b, -branch           optimize branch lengths as well 
  -o, -output           output parameter file 
  ```
  
  The simulation program only generates global jumps. To convert global jumps to
  local paths, users need to run below command:
  ```
  ./bin/global_jumps_to_paths [OPTIONS] <treefile> <statefile> <global jumps file> <output local paths>
  ```
  Available options include:
  ```
  -v, -verbose  print more run info 
  ```

### A simple test run
```
cd test
./run.sh
```
The above command will run a test including three steps:
1. Simulate one evolution path of 10-site long sequence.
2. Fix the starting and ending states, collect 5000 simulated paths and print summary statistics.
3. Run a MCMC sampler to get 5000 paths and print summary statistics. A report will be generated in `figures/N10`.

Below is the complete usage of `run.sh`:
```
Usage: run.sh
          [-n MCMC-EM iterations (5000) ] [-s sites (10)]
          [-P proposal 0:poisson 1:direct 2:forward 3:unif (0)]
          [-T sample full tree (false)] [-E estimate parameters (false)]
          [-L MCMC burn-in (0)] [-B MCMC batch (10)]
          [-R fix root (false)] [-r root sequence (null)]
          [-f output file prefix (test)]
          [-p true parameter file (input/test.param)]
          [-t true tree (input/test.nwk)]
          [-i initial parameter file (input/test.param)]
          [-j initial tree (input/test.nwk)]
          [-k initial (local) path file (null)]
```


Contacts
========================

Andrew D. Smith
andrewds@usc.edu

Xiaojing Ji
xiaojinj@usc.edu


Copyright and License Information
=================================

Copyright (C) 2018-2020
University of Southern California,
Andrew D. Smith
  
Current Authors:  Andrew D. Smith, Xiaojing Ji and Jenny Qu
  
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
  
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
  
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
