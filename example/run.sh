prefix=test
../src/prog/epievo_sim ${prefix}.param tree.nwk -o ${prefix}.outmeth -n 100000 -p ${prefix}.global_jumps -v 
../src/prog/global_jumps_to_paths tree.nwk ${prefix}.outmeth ${prefix}.global_jumps ${prefix}.paths
Rscript plot.R ${prefix} 2 0 5000 



