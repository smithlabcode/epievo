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


Usage
========================
### Simulate epigenome evolution
`epievo_sim` can be used to simulate epigenome evolution from given evolution parameters,
tree topology, or the total time duration of a single branch. 
#### Example
```
cd test
../bin/epievo_sim -v -n 1000 -o tree.states -p tree.global_jumps -t tree.nwk tree.param
../bin/global_jumps_to_paths tree.nwk tree.states tree.global_jumps tree.local_paths
```
Two output files will be generated from `epievo_sim`:
`tree.global_jumps` contains mutation information ordered by position
and time, and `tree.states` contains epigenomic states at each position and each node.
The second command will convert `tree.global_jumps` to `tree.local_paths`,
which contains mutation information organized at each single site.
Local paths are used for model inference.

Below is the usage of `epievo_sim`:
```
Usage: epievo_sim [OPTIONS] <params-file>

Options:
  -o, -output    name of output file for epigenomic states(default: stdout) 
  -n, -n-sites   length of sequence to simulate (default: 100) 
  -p, -paths     name of output file for evolution paths as sorted jump times 
                 (default: stdout) 
  -s, -seed      rng seed 
  -r, -root      root states file 
  -t, -tree      Newick format tree file 
  -T, -evo-time  evolutionary time 
  -l, -leaf      write only leaf states (default: false) 
  -S, -scale     scale the branch lengths (default: true 
  -R, -rates     use triplet transition rates (default: false) 
  -v, -verbose   print more run info 

Help options:
  -?, -help      print this help message 
      -about     print about message 
```

### Estimate model parameters from complete history
`bin/epievo_est_complete` is used to obtain maximum-likelihood estimates
of model parameters (and branch lengths if specified),
from provided local paths of complete evolution history.
#### Example
```
cd test
../bin/epievo_est_complete -v -p tree.param -t tree.nwk -o tree.param.updated tree.local_paths
```
The output file `tree.param.updated` contains estimated model parameters
(and branch lengths if specified).
Below is the usage of `epievo_est_complete`:
```
Usage: epievo_est_complete [OPTIONS] <path-file>

Options:
  -p, -param    initial parameter file 
  -t, -tree     initial tree file in newick format 
  -v, -verbose  print more run info 
  -b, -branch   optimize branch lengths as well 
  -o, -output   output parameter file 

Help options:
  -?, -help     print this help message 
      -about    print about message 
  ```

### Estimate model parameters and histories from leaf data
Program `epievo_est_params_histories` runs a MCMC-EM algorithm to estimate model
parameters and sample evolution histories simultaneously, which requires
initial parameters, local paths to be provided.
If you only have observed data (epigenomic states at leaf species), you can
use `epievo_initialization` to generate starting parameters and paths
through heuristics and site-independent-model-based methods.

#### Example
```
../bin/epievo_initialization -p tree.param.init -o tree.local_paths.init tree.nwk observed.states
../bin/epievo_est_params_histories -v -o tree.local_path.est -p tree.local_path.est \
  tree.param.init tree.nwk tree.local_paths.init
```

Below is the usage of `epievo_est_params_histories`:
```
Usage: epievo_est_params_histories [OPTIONS] <param> <treefile> <path_file>

Options:
  -i, -iteration  number of MCMC-EM iteration (default: 10) 
  -B, -batch      number of MCMC iteration (default: 10) 
  -L, -burnin     MCMC burn-in length (default: 10) 
  -s, -seed       rng seed 
  -o, -outfile    output file of local paths 
  -p, -outparam   output file of parameters 
  -b, -branch     optimize branch lengths as well 
  -v, -verbose    print more run info 

Help options:
  -?, -help       print this help message 
      -about      print about message 
```

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
Local path is another way to organize mutation events.
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
