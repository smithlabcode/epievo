num=5000
proposal=2
burnin=1000
batch=10

estParam=0
sampleTree=0
fixRoot=0

paramFile=""
initPathFile=""
initTreeFile=""
outPrefix="out"
outDir=$PWD

print_usage() {
  printf "Usage: $(basename $0) [-n MCMC-EM iterations] [-P proposal]
          [-L MCMC burn-in length] [-B MCMC batch size]
          [-E estimate parameters] [-T sample full tree] [-R fix root]
          [-p parameter file] [-i initial path file] [-t initial tree (.nwk)]
          [-f output file prefix] [-o output directory]\n"
}

while getopts 'n:P:L:B:ETRp:i:t:f:o:h' flag; do
  case "${flag}" in
    n) num="${OPTARG}" ;;
    P) proposal="${OPTARG}" ;;
    L) burnin="${OPTARG}" ;;
    B) batch="${OPTARG}" ;;
    E) estParam=1 ;;
    T) sampleTree=1 ;;
    R) fixRoot=1 ;;
    p) paramFile="${OPTARG}" ;;
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


outPathFile=$outDir/$outPrefix.update.path_local
statsFile=$outDir/$outPrefix.stats
traceFile=$outDir/$outPrefix.trace

#------------------------------------------------------------------------------
CMD="mcmc_test -P $proposal -B $batch -i $num -L $burnin
-o $outPathFile -S $statsFile -t $traceFile
$paramFile $initTreeFile $initPathFile"

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
