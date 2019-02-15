num=50
sites=10000

outPrefix=test

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

paramFile=$inDir/$outPrefix.param
treeFile=$inDir/$outPrefix.nwk

################################################################################

print_usage() {
  printf "Usage: $(basename $0)
          [-n number of samples (50) ] [-s sites (10000)]
          [-f output file prefix (test)]
          [-p parameter file (input/test.param)]
          [-t tree file (input/test.nwk)]\n"
}

while getopts 'n:s:f:p:t:h' flag; do
  case "${flag}" in
    n) num="${OPTARG}" ;;
    s) sites="${OPTARG}" ;;
    f) outPrefix="${OPTARG}" ;;
    p) paramFile="${OPTARG}" ;;
    t) treeFile="${OPTARG}" ;;
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
  epievo_sim $paramFile -t $treeFile -n $sites \
    -o $outDir/$outPrefix.outmeth -p $outDir/$outPrefix.global_jumps
  global_jumps_to_paths $treeFile $outDir/$outPrefix.outmeth \
    $outDir/$outPrefix.global_jumps $outDir/$outPrefix.path_local
  epievo_est_complete -S -p $paramFile -t $treeFile \
    -o $outDir/$outPrefix.param.update $outDir/$outPrefix.path_local
 
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

Rscript ../rscripts/compare_params.R -t $paramFile -e $outDir/$outPrefix \
  -o $figDir/${outPrefix}_N$sites.pdf
