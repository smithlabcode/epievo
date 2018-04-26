prefix=test
# simulate evolution
../src/prog/epievo_sim -s 1 ${prefix}.param tree.nwk -o ${prefix}.outmeth -n 100000 -p ${prefix}.global_jumps -v 
# convert format from global path to local paths
../src/prog/global_jumps_to_paths tree.nwk ${prefix}.outmeth ${prefix}.global_jumps ${prefix}.paths
# plot local paths of branch in range
Rscript plot.R ${prefix} 2 0 5000 
# estimate model parameter from local paths
../src/prog/epievo_est_complete -p ${prefix}.param -t tree.nwk -v -b ${prefix}.paths



