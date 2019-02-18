num=100
burnin=10
batch=10

estParam=0
sampleTree=0
fixRoot=0

initParamFile=input/test.param
initPathFile=output/sim_results/test.path_local
initTreeFile=input/test.nwk
outPrefix=test
outDir=output/mcmc/poisson

print_usage() {
  printf "Usage: $(basename $0)
          [-n MCMC-EM iterations (10)]
          [-L MCMC burn-in (10)] [-B MCMC batch (10)]
          [-E estimate parameters (false)]
          [-T sample full tree (false)] [-R fix root (false)]
          [-p initial parameter file (input/test.param)]
          [-i initial path file (output/sim_results/test.path_local)]
          [-t initial tree (input/test.nwk)]
          [-f output file prefix (test)]
          [-o output directory (output/mcmc/poisson)]\n"
}

while getopts 'n:L:B:ETRp:i:t:f:o:h' flag; do
  case "${flag}" in
    n) num="${OPTARG}" ;;
    L) burnin="${OPTARG}" ;;
    B) batch="${OPTARG}" ;;
    E) estParam=1 ;;
    T) sampleTree=1 ;;
    R) fixRoot=1 ;;
    p) initParamFile="${OPTARG}" ;;
    i) initPathFile="${OPTARG}" ;;
    t) initTreeFile="${OPTARG}" ;;
    f) outPrefix="${OPTARG}" ;;
    o) outDir="${OPTARG}" ;;
    h) print_usage
       exit 0 ;;
    *) print_usage
       exit 1 ;;
  esac
done

mkdir -p $outDir
outPathFile=$outDir/$outPrefix.update.path_local
statsFile=$outDir/$outPrefix.stats
traceFile=$outDir/$outPrefix.trace

#------------------------------------------------------------------------------
CMD="mcmc_test -B $batch -i $num -L $burnin
-o $outPathFile -S $statsFile -t $traceFile
$initParamFile $initTreeFile $initPathFile"

if [ "$estParam" -eq 1 ]; then
  CMD="${CMD} -m"
fi

if [ "$fixRoot" -eq 1 ]; then
  CMD="${CMD} -R"
fi

if [ "$sampleTree" -eq 0 ]; then
  CMD="${CMD} -b"
fi

echo $CMD
eval $CMD
