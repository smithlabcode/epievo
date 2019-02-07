proposalTexts=( 'poisson' 'direct' 'forward' 'unif' )
proposalToTest=( 0 1 2 3 )

num=5000
sites=10
burnin=0
batch=10

proposal=0
sampleTree=0
estParam=0
fixRoot=0

rootseq=""
outPrefix=test

################################################################################
####### ##    ###   ##    ###    ##   ###    #    ##    ###   ###    ##
####### # #    #    # #   #     #      #    # #   # #    #    #     #
####### # #    #    ##    ##    #      #    # #   ##     #    ##     #
####### # #    #    # #   #     #      #    # #   # #    #    #       #
####### ##    ###   # #   ###    ##    #     #    # #   ###   ###   ##

inDir=input
outDir=output
simDir=$outDir/sim_results
fwDir=$outDir/forward
mcmcRootDir=$outDir/mcmc
figDir=figures

mkdir -p $simDir
mkdir -p $fwDir
mkdir -p $mcmcRootDir
mkdir -p $figDir
################################################################################

trueParamFile=$inDir/$outPrefix.param
trueTreeFile=$inDir/$outPrefix.nwk
initPathFile=$simDir/$outPrefix.path_local
initParamFile=$trueParamFile
initTreeFile=$trueTreeFile

################################################################################


print_usage() {
  printf "Usage: $(basename $0) [-n MCMC-EM iterations] [-s sites]
          [-P proposal (default: poisson)] [-T sample full tree]
          [-E estimate parameters]
          [-L burnin] [-B batch ] [-R fix root] [-r root sequence]
          [-f output file prefix] [-p true parameter file] [-t true tree (.nwk)]
          [-i initial parameter file] [-j initial tree (.nwk)]
          [-k initial (local) path file]\n"
}

while getopts 'n:s:P:TEL:B:Rr:f:p:t:i:j:k:h' flag; do
  case "${flag}" in
    n) num="${OPTARG}" ;;
    s) sites="${OPTARG}" ;;
    P) proposal="${OPTARG}" ;;
    T) sampleTree=1 ;;
    E) estParam=1
       batch=10;;
    L) burnin="${OPTARG}" ;;
    B) batch="${OPTARG}" ;;
    R) fixRoot=1 ;;
    r) rootseq="${OPTARG}" ;;
    f) outPrefix="${OPTARG}" ;;
    p) trueParamFile="${OPTARG}"
       initParamFile=$trueParamFile;;
    t) trueTreeFile="${OPTARG}"
       initTreeFile=$trueTreeFile;;
    i) initParamFile="${OPTARG}" ;;
    j) initTreeFile="${OPTARG}" ;;
    k) initPathFile="${OPTARG}" ;;
    h) print_usage
       exit 0 ;;
    *) print_usage
       exit 1 ;;
  esac
done

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################



#------------------------------------------------------------------------------
echo "------------------------------------------------------------------------"
echo Number of samples: $num
echo Number of sites: $sites
echo Sample full tree: $sampleTree

########### STEP 1 ############################################################
echo "------------------------------------------------------------------------"
echo 1. Generate initial paths...

simCMD="./sim.sh -n 1 -s $sites -f $outPrefix -o $simDir
-p $trueParamFile -t $trueTreeFile"
if [ "$rootseq" != "" ]; then
  simCMD="${simCMD} -r $rootseq"
fi
eval $simCMD

########### STEP 2 ############################################################
if [ "$estParam" -eq 0 ]; then
  echo "------------------------------------------------------------------------"
  echo 2. End-conditioned forward sampling...
  forwardCMD="./forward.sh -n $num -s $sites 
-p $trueParamFile -t $trueTreeFile -i $initPathFile -f $outPrefix -o $fwDir"
  if [ "$sampleTree" -eq 1 ]; then
    forwardCMD="${forwardCMD} -T"
  fi
  if [ "$fixRoot" -eq 1 ]; then
    forwardCMD="${forwardCMD} -R"
  fi
  eval $forwardCMD
fi

############ STEP 3 ###########################################################
echo "------------------------------------------------------------------------"
echo 3. MCMC...
echo Burn-in: $burnin
echo Batch: $batch
echo Output directory: $mcmcRootDir

if [ "$proposal" -lt ${#proposalTexts[@]} ]; then
  proposalToTest=( ${proposalToTest[$proposal]} )
fi

for (( i=0; i<${#proposalToTest[@]}; i++ )); do
  proposal=${proposalToTest[$i]}
  proposalLabel=${proposalTexts[$proposal]}
  echo Proposal: $proposalLabel
   
  mcmcDir=$mcmcRootDir/$proposalLabel
  reportDir=$figDir/N$sites
  mkdir -p $mcmcDir
  mkdir -p $reportDir

  #--------------------------------------------------------------------------
  mcmcCMD="./mcmc.sh -n $num -P $proposal -L $burnin -B $batch
-p $initParamFile -i $initPathFile -t $initTreeFile
-f $outPrefix -o $mcmcDir"

  if [ "$sampleTree" -eq 1 ]; then
    mcmcCMD="${mcmcCMD} -T"
  fi
  if [ "$estParam" -eq 1 ]; then
    mcmcCMD="${mcmcCMD} -E"
  fi
  if [ "$fixRoot" -eq 1 ]; then
    mcmcCMD="${mcmcCMD} -R"
  fi
  eval $mcmcCMD


  if [ "$estParam" -eq 0 ]; then
    Rscript ../rscripts/sufficient_stats_distr.R \
      --true $fwDir/$outPrefix \
      --est $mcmcDir/$outPrefix \
      --output $reportDir/${outPrefix}_${proposalLabel}.pdf \
      ${proposalToTest[@]}
  fi
done

if [ "$estParam" -eq 1 ]; then
  Rscript ../rscripts/sufficient_stats_conv.R \
    --prefix $outPrefix --input $mcmcRootDir \
    --output $reportDir/${outPrefix}_est.pdf \
    --param $trueParamFile \
    ${proposalToTest[@]}
fi

echo "------------------------------------------------------------------------"
