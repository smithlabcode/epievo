proposalTexts=( 'direct' 'unif' 'poisson' 'forward' )
proposalToTest=( 0 1 2 3 )
#proposalToTest=( 2 )

num=1000
sites=5
rootseq=""
sampleTree=0
testConv=0
fixRoot=0
prefix="test"
proposal=${#proposalTexts[@]}
burning=1000
batch=1

print_usage() {
  printf "Usage: $(basename $0) [-n number] [-s sites] [-r root sequence]
          [-P proposal] [-f file prefix] [-T sample full tree]
          [-C track convergence]
          [-L burning] [-B batch ] [-R fix root]\n"
}

while getopts 'n:s:r:P:f:TCL:B:Rh' flag; do
  case "${flag}" in
    n) num="${OPTARG}" ;;
    s) sites="${OPTARG}" ;;
    r) rootseq="${OPTARG}" ;;
    P) proposal="${OPTARG}" ;;
    f) prefix="${OPTARG}" ;;
    T) sampleTree=1 ;;
    C) testConv=1
       burning=0
       batch=10;;
    L) burning="${OPTARG}" ;;
    B) batch="${OPTARG}" ;;
    R) fixRoot=1 ;;
    h) print_usage
       exit 0 ;;
    *) print_usage
       exit 1 ;;
  esac
done


dpar=input
dout=output
dsim=$dout/sim_results
dfw=$dout/forward
dlog=log

#------------------------------------------------------------------------------
echo "------------------------------------------------------------------------"
echo Number of samples: $num
echo Number of sites: $sites
echo File prefix: $prefix
echo Sample full tree: $sampleTree

########### STEP 1 ############################################################
echo "------------------------------------------------------------------------"
echo 1. Generate initial paths...
simCMD="./sim.sh -n 1 -s $sites -f $prefix -p $dpar -o $dsim"
if [ "$rootseq" != "" ]; then
  simCMD="${simCMD} -r $rootseq"
fi
eval $simCMD

########### STEP 2 ############################################################
if [ "$testConv" -eq 0 ]; then
  echo 2. End-conditioned forward sampling...
  forwardCMD="./forward.sh -n $num -s $sites -f $prefix -p $dpar -i $dsim -o $dfw"
  if [ "$fixRoot" -eq 1 ]; then
    forwardCMD="${forwardCMD} -R"

  fi
  eval $forwardCMD
fi

########### STEP 3 ############################################################
echo "------------------------------------------------------------------------"
echo 3. MCMC...
echo Burning: $burning
echo Batch: $batch
if [ "$proposal" -lt ${#proposalTexts[@]} ]; then
  proposalToTest=( ${proposalToTest[$proposal]} )
fi

for (( i=0; i<${#proposalToTest[@]}; i++ )); do
  proposal=${proposalToTest[$i]}
  proposalLabel=${proposalTexts[$proposal]}
  echo Proposal: $proposalLabel
 
  dmcmc=$dout/mcmc/$proposalLabel
  dreport=figures/N${sites}
  mkdir -p $dmcmc
  mkdir -p $dreport

  #--------------------------------------------------------------------------
  mcmcCMD="./mcmc.sh -n $num -s $sites -f $prefix -p $dpar -i $dsim -o $dmcmc
-P $proposal -L $burning -B $batch"
  if [ "$sampleTree" -eq 1 ]; then
    mcmcCMD="${mcmcCMD} -T"
  fi
  if [ "$testConv" -eq 1 ]; then
    mcmcCMD="${mcmcCMD} -E"
  fi
  if [ "$fixRoot" -eq 1 ]; then
    mcmcCMD="${mcmcCMD} -R"
  fi
  mcmcCMD="${mcmcCMD} > log/log_${prefix}_${proposal}_N${sites}"
  eval $mcmcCMD


  if [ "$testConv" -eq 0 ]; then
    sufficient_stats_distr.R \
      --true $dfw/${prefix}_N${sites} \
      --est $dmcmc/${prefix}_N${sites} \
      --output $dreport/${prefix}_N${sites}_${proposalLabel}.pdf \
      ${proposalToTest[@]}
  fi
done

if [ "$testConv" -eq 1 ]; then
  sufficient_stats_conv.R \
    --prefix ${prefix}_N${sites} --input $dout/mcmc \
    --output $dreport/${prefix}_N${sites}_all.pdf \
    --param $dpar/$prefix.param \
    ${proposalToTest[@]}
fi

