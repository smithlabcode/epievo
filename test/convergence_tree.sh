proposalTexts=( 'direct' 'unif' 'poisson' 'forward' )
proposalToTest=( 2 )
num=${1:-1000}
sites=${2:-5}
batch=${3:-1}
proposal=${4:-${#proposalTexts[@]}}


cd ~/data/evo/mcmc

prefix=tree
dpar=input
dout=output
dsim=$dout/sim_results
dfw=$dout/forward
dlog=log
dreport=figures/N${sites}
mkdir -p $dreport

echo Number of samples: $num
echo Number of sites: $sites

echo 1. Generate initial paths...
./sim.sh 1 $sites $prefix $dpar $dsim

echo 2. MCMC...
if [ "$proposal" -lt ${#proposalTexts[@]} ]; then
  proposalTexts=( ${proposalTexts[$i]} )
fi

for (( i=0; i<${#proposalToTest[@]}; i++ )); do
  proposal=${proposalToTest[$i]}
  proposalLabel=${proposalTexts[$proposal]}
  echo Proposal: $proposalLabel
  
  dmcmc=$dout/mcmc/$proposalLabel
  mkdir -p $dmcmc

  #--------------------------------------------------------------------------

  ./mcmc_tree.sh $num $sites $prefix $dpar $dsim $dmcmc $proposal 1 0 $batch\
    > log/log_N${sites}_$proposalLabel

done
 
~/data/evo/rscripts/test_mcmc/sufficient_stats_conv.R \
        --prefix ${prefix}_N${sites} --input $dout/mcmc \
        --output $dreport/${prefix}_N${sites}_all.pdf \
        --param $dpar/$prefix.param \
        ${proposalToTest[@]}
