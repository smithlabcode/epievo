num=$1
sites=$2
prefix=$3
dpar=$4
dsim=$5
dmcmc=$6
proposal=$7
estparam=${8:-0}
burning=${9:-1000} 
batch=${10:-1} 


dsrc=~/codes/EvoSim/src                                                         
outprefix=${prefix}_N${sites}

param=$dpar/$prefix.param
tree=$dpar/$prefix.nwk
inpaths=$dsim/$outprefix.path_local
outpaths=$dmcmc/$outprefix.update.path_local
stats=$dmcmc/$outprefix.stats
trace=$dmcmc/$outprefix.trace

#------------------------------------------------------------------------------

if [ "$estparam" -gt 0 ]; then
$dsrc/harnesses/metropolis_unif_test -P $proposal -B $batch -r $num -L $burning \
        -S $stats -t $trace \
        -o $outpaths -m \
        $param $tree $inpaths
else
$dsrc/harnesses/metropolis_unif_test -P $proposal -B $batch -r $num -L $burning \
        -S $stats -t $trace \
        -o $outpaths \
        $param $tree $inpaths
fi
