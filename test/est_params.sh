num=50
sites=10000

estRatesTree=0
outPrefix=tree

################################################################################
####### ##    ###   ##    ###    ##   ###    #    ##    ###   ###    ##
####### # #    #    # #   #     #      #    # #   # #    #    #     #
####### # #    #    ##    ##    #      #    # #   ##     #    ##     #
####### # #    #    # #   #     #      #    # #   # #    #    #       #
####### ##    ###   # #   ###    ##    #     #    # #   ###   ###   ##

inDir=input
outDir=output/mle
figDir=figures/mle

mkdir -p $outDir
mkdir -p $figDir

################################################################################

trueParamFile=$inDir/$outPrefix.param
trueTreeFile=$inDir/$outPrefix.nwk
initParamFile=$trueParamFile
initTreeFile=$trueTreeFile


################################################################################

print_usage() {
  printf "Usage: $(basename $0)
          [-n number of samples (50) ] [-s sites (10000)]
          [-f output file prefix (tree)]
          [-T estimate both rates and branch lengths]
          [-p parameter file (input/tree.param)]
          [-t tree file (input/tree.nwk)]
          [-i initial parameter file (input/tree.param)]
          [-j initial tree (input/tree.nwk)]\n"
}

while getopts 'n:s:f:Tp:t:i:j:h' flag; do
  case "${flag}" in
    n) num="${OPTARG}" ;;
    s) sites="${OPTARG}" ;;
    f) outPrefix="${OPTARG}" ;;
    T) estRatesTree=1 ;; 
    p) trueParamFile="${OPTARG}"
       initParamFile=$trueParamFile ;;
    t) trueTreeFile="${OPTARG}"
       initTreeFile=$trueTreeFile ;;
    i) initParamFile="${OPTARG}" ;;
    j) initTreeFile="${OPTARG}" ;;
    h) print_usage
       exit 0 ;;
    *) print_usage
       exit 1 ;;
  esac
done

# estimates of parameters
> $outDir/${outPrefix}_st0.stats
> $outDir/${outPrefix}_st1.stats
> $outDir/${outPrefix}_bl0.stats
> $outDir/${outPrefix}_bl1.stats
> $outDir/${outPrefix}_rt0.stats
> $outDir/${outPrefix}_rt1.stats

for (( i=0; i<$num; i++ )); do
  epievo_sim $trueParamFile -t $trueTreeFile -n $sites \
    -o $outDir/$outPrefix.outmeth -p $outDir/$outPrefix.global_jumps
  global_jumps_to_paths $trueTreeFile $outDir/$outPrefix.outmeth \
    $outDir/$outPrefix.global_jumps $outDir/$outPrefix.path_local
  mleCMD="epievo_est_complete -S -p $initParamFile -t $initTreeFile -v
-o $outDir/$outPrefix.param.update $outDir/$outPrefix.path_local"
  if [ "$estRatesTree" -eq 1 ]; then
      mleCMD="${mleCMD} -b"
  fi
  eval $mleCMD
 
  updatedParamFile=$outDir/$outPrefix.param.update
  text=$(<$updatedParamFile)
  st0=`echo $text | cut -d' ' -f2`
  st1=`echo $text | cut -d' ' -f3`
  bl0=`echo $text | cut -d' ' -f5`
  bl1=`echo $text | cut -d' ' -f6`
  rt0=`echo $text | cut -d' ' -f8`
  rt1=`echo $text | cut -d' ' -f9`

  echo $st0 >> $outDir/${outPrefix}_st0.stats
  echo $st1 >> $outDir/${outPrefix}_st1.stats
  echo $bl0 >> $outDir/${outPrefix}_bl0.stats
  echo $bl1 >> $outDir/${outPrefix}_bl1.stats
  echo $rt0 >> $outDir/${outPrefix}_rt0.stats
  echo $rt1 >> $outDir/${outPrefix}_rt1.stats

done

Rscript ../rscripts/compare_params.R -t $trueParamFile -e $outDir/$outPrefix \
  -o $figDir/${outPrefix}_N$sites.pdf
