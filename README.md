Epievo includes tools for simulation and inference of epigenome evolution.
The epigenome evolution is modeled by a context-dependent continuous-time
Markov process, and estimated based
on a Monte Carlo Expectation-Maximization algorithm.


Building and Installing 
=======================

First, download `epievo` in recursive mode,
```
git clone --recursive https://github.com/smithlabcode/epievo.git
```
then install `epievo` by running
```
cd epievo
make install
```

Usage
========================
### Simulating epigenome evolution
`epievo_sim` takes [model parameters](#model-parameters) to
simulate epigenomic states of each sepecies
and internal jumps, based on a context-dependent continuous-time
Markov model. It will create two outputs:
[epigenomic states](#epigenomic-states) at each genomic position and each
node,
and [global jumps](#global-jumps) ordered by node label, genomic
position and evolutionary time (file name could be provided after argument `-p`).

```
epievo_sim [options] <parameter file> <outfile>
```
Users can use `epievo_sim` to simulate multiple species (using argument `-t` to
provide a [phylogenetic tree](#tree-format)), or one single branch by using
argument `-T` to set total evolution time.
Epigenomic states of root species can be fixed by providing a states file
after `-r`. By default, `epievo_sim` will normalize input parameters following
"one-mutation-per-unit-time-per-site" rule. To use un-normalized
parameters, users can use flat `-S`.

### Converting global jumps to local paths
The [global jumps](#global-jumps) data structure allows fast forward simulation,
and ([local paths](#local-paths) data structure is more efficient and used in
backward history inference.
Global-jump files can be converted to local paths through program `global_jumps_to_paths`:
```
global_jumps_to_paths [options] <statefile> <jumpfile> <outfile>
```
Users can pass [Phylogenetic tree](#tree-format) in newick format after argument `-t`,
or set the total evolution time of a single branch after argument `-T`.


### Estimating model parameters from complete history
`epievo_est_complete` is used to estimate [model parameters](#model-parameters)
when the complete information of evolution process is given
([evolution paths](#local-paths) and [tree](tree-format)) are known.
The maximum-likelihood estimates are obtained based on gradient-ascent approach.
Initial values of model parameters are required.
By default, `epievo_est_complete` will not update tree branches. To
learn model parameters and branches together, flag `-b` needs to be specified.
```
epievo_est_complete [options] <parameter file> <tree file> <local paths>
```

### Estimating model parameters and histories from leaf data
In practice, epigenomic states are only observed in leaf species.
Programs `epievo_initialization` and `epievo_est_params_histories` allow
users to estimate [evolution paths](#local-paths) and
[model parameters](#model-parameters) simultaneously, given
the leaf data and a starting tree (e.g. setting all branches to ones). 

`epievo_initialization` is used to generate initial evolution histories
and model parameters through heuristics and site-independent-model-based
methods.
```
epievo_initialization [options] <tree file> <states file>

```

Program `epievo_est_params_histories` runs a MCEM algorithm to estimate model
parameters and sample evolution histories iteratively, which requires
initial parameters, [local paths](#local-paths) to be provided.
By default, only model parameters will be estimated and printed to output
file (specified by `-p`).
To estimate branch lengths simultanesously, 
users need to pass the `-b` flag.
If only one branch is included in the data, users should pass the `-T` flag.
Other training parameters include MCEM total iterations (`-i`),
MCMC sample size (`-B`) and MCMC burn-in length (`-L`).
```
epievo_est_params_histories [options] <parameter file> <tree file/time> <local paths>
```

### Inferring histories between two given state-sequences
Program `epievo_sim_pairwise` runs a MCMC algorithm to infer epigenomic evolution
between two given state-sequences. The output will be [local paths](#local-paths)
between ending sequences.
```
epievo_sim_pairwise [OPTIONS] <parameter file> <states file>
```
If only one branch is included in the data, users should pass the `-T` flag.
MCMC burn-in length can be specified after argument `-L`.


Running the tests
========================

The command below will generate the complete evolution information from
a phylogenetic tree in `tree.nwk`, and model parameters in `test.param`.
```
cd test
../bin/epievo_sim -v -n 1000 -o test.states -p test.global_jumps -t tree.nwk test.param
```
Two output files will be generated from `epievo_sim`.
`test.global_jumps` contains mutation information ordered by position
and time, and `test.states` contains epigenomic states at each position and each node.

To run inference programs, the global jumps should be converted to
local paths first, by running:
```
../bin/global_jumps_to_paths -v tree.nwk test.states test.global_jumps test.local_paths
```

The command below will estimate model parameters (saved in `test.param.updated`) from
tree file `tree.nwk`, local paths `test.local_paths`, given starting parameters
`test.param`.
```
../bin/epievo_est_complete -v -o test.param.updated test.param tree.nwk test.local_paths
```

Now, we can try to initialize the inference procedure from a tree `tree.nwk` and
leaf data `observed.states`. Initial estimates of parameters and evolution histories
will be saved in `test.param.init` and `test.local_paths.init` respectively.
```
../bin/epievo_initialization -v -p test.param.init -o test.local_paths.init tree.nwk observed.states
```
Then, the command below will run a MCEM procedure to estimate model parameters and
evolution histories, which will be printed in `test.local_path.est`
and `test.local_path.est` respectively.
```
../bin/epievo_est_params_histories -v -o test.local_path.est -p test.local_path.est \
  test.param.init tree.nwk test.local_paths.init
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

### Local paths
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
