num=$1
sites=$2
prefix=$3
dpar=$4
dsim=$5

dsrc=~/codes/EvoSim/src
outprefix=${prefix}_N${sites}

param=$dpar/$prefix.param
tree=$dpar/$prefix.nwk
outmeth=$dsim/$outprefix.outmeth
jumps=$dsim/$outprefix.global_jumps
paths=$dsim/$outprefix.path_local
stats=$dsim/$outprefix.stats

#------------------------------------------------------------------------------

$dsrc/harnesses/stationary_stats_est $param -t $tree \
        -n $sites -m $num \
        -o $outmeth -p $jumps -S $stats

$dsrc/prog/global_jumps_to_paths $tree $outmeth \
        $jumps $paths
