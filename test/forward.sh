num=1000
sites=5

paramFile=""
treeFile=""
initPathFile=""
outPrefix="out"
outDir=$PWD

sampleTree=0
fixRoot=0

print_usage() {
  printf "Usage: $(basename $0) [-n number] [-s sites]
          [-p parameter file] [-t tree (.nwk)] [-i input path file]
          [-f file prefix] [-o output directory]
          [-T sample full tree] [-R fix root]\n"
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

