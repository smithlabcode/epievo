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
To run the full pipeline,
`Python2.x` and `R` need to be installed and added to environment variables.
You will also need to have R package `optparse`,
`ggplot2`
and Python libraries `argparse` installed.

Usage
========================
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
          [-n MCMC-EM iterations (10) ] [-s sites (10)]
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
