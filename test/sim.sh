num=1
sites=5

paramFile=""
treeFile=""
outPrefix="out"
outDir=$PWD

rootseq=""

print_usage() {
  printf "Usage: $(basename $0) [-n number] [-s sites] [-r root sequence]
          [-p parameter file] [-t tree (.nwk)]
          [-f output file prefix] [-o output directory]\n"
}

while getopts 'n:s:r:p:t:f:o:h' flag; do
  case "${flag}" in
    n) num="${OPTARG}" ;;
    s) sites="${OPTARG}" ;;
    r) rootseq="${OPTARG}" ;;
    p) paramFile="${OPTARG}" ;;
    t) treeFile="${OPTARG}" ;;
    f) outPrefix="${OPTARG}" ;;
    o) outDir="${OPTARG}" ;;
    h) print_usage
       exit 0 ;;
    *) print_usage
       exit 1 ;;
  esac
done

inDir=$(dirname $paramFile)
rootFile=$inDir/$rootseq.rootseq

statesFile=$outDir/$outPrefix.outmeth
globalPathFile=$outDir/$outPrefix.global_jumps
pathFile=$outDir/$outPrefix.path_local

#------------------------------------------------------------------------------
if [ "$rootseq" = "" ]; then
  epievo_sim $paramFile -t $treeFile \
    -v -n $sites -o $statesFile -p $globalPathFile
else
  python ../pyscripts/write_rootseq_to_file.py $rootseq $rootFile

  epievo_sim $paramFile -t $treeFile -r $rootFile \
    -v -n $sites -o $statesFile -p $globalPathFile
fi

global_jumps_to_paths $treeFile $statesFile $globalPathFile $pathFile
