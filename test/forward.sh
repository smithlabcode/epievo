num=5000
sites=10

paramFile=input/test.param
treeFile=input/test.nwk
initPathFile=output/sim_results/test.path_local
outPrefix=test
outDir=output/forward

sampleTree=0
fixRoot=0

print_usage() {
  printf "Usage: $(basename $0)
          [-n number (5000)] [-s sites (10)]
          [-p parameter file (input/test.param)]
          [-t tree (input/test.nwk)]
          [-i input path file (output/sim_results/test.path_local)]
          [-f file prefix (test)] [-o output directory (output/forward)]
          [-T sample full tree (false)] [-R fix root (false)]\n"
}

while getopts 'n:s:p:t:i:f:o:TRh' flag; do
  case "${flag}" in
    n) num="${OPTARG}" ;;
    s) sites="${OPTARG}" ;;
    p) paramFile="${OPTARG}" ;;
    t) treeFile="${OPTARG}" ;;
    i) initPathFile="${OPTARG}" ;;
    f) outPrefix="${OPTARG}" ;;
    o) outDir="${OPTARG}" ;;
    T) sampleTree=1 ;;
    R) fixRoot=1 ;;
    h) print_usage
       exit 0 ;;
    *) print_usage
       exit 1 ;;
  esac
done

mkdir -p $outDir
pathFile=$outDir/$outPrefix.path_local
statsFile=$outDir/$outPrefix.stats

#-------------------------------------------------------------------------------
CMD="forward_sim_stats -n $num -S $statsFile -o $pathFile
$paramFile $treeFile $initPathFile"

if [ "$fixRoot" -eq 1 ]; then
  CMD="${CMD} -R"
fi

if [ "$sampleTree" -eq 0 ]; then
  CMD="${CMD} -b"
fi

eval $CMD

