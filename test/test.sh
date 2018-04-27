prefix=test
../src/prog/epievo_sim  -s 1 ${prefix}.param tree.nwk -o ${prefix}.outmeth -n 10 -p ${prefix}.global_jumps -v 
../src/prog/global_jumps_to_paths tree.nwk ${prefix}.outmeth ${prefix}.global_jumps ${prefix}.path_local
../src/prog/test_sampling -s 1 -p ${prefix}.param -t tree.nwk -r 5 -o ${prefix}.update.path_local -S ${prefix}.update.outmeth ${prefix}.path_local
../src/prog/epievo_est_complete -p ${prefix}.param -t tree.nwk  -v -b ${prefix}.path_local
../src/prog/epievo_est_complete -p ${prefix}.param -t tree.nwk  -v -b ${prefix}.update.path_local

for r in `echo 1 2 5 10`; do
echo $r;
../src/prog/test_sampling -s 1 -p ${prefix}.param -t tree.nwk -r ${r} -o ${prefix}.update_${r}.path_local -S ${prefix}.update.outmeth ${prefix}.path_local
done
