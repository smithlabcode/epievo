sites=10

paramFile=input/test.param
treeFile=input/test.nwk
outPrefix=test
outDir=output/sim_results

rootseq=""

print_usage() {
  printf "Usage: $(basename $0)
          [-s sites (10)] [-r root sequence]
          [-p parameter file (input/test.param)]
          [-t tree (input/test.nwk)]
          [-f output file prefix (test)]
          [-o output directory (output/sim_results)]\n"
}

while getopts 's:r:p:t:f:o:h' flag; do
  case "${flag}" in
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

mkdir -p $outDir
inDir=$(dirname $paramFile)
rootFile=$inDir/$rootseq.rootseq

statesFile=$outDir/$outPrefix.outmeth
globalPathFile=$outDir/$outPrefix.global_jumps
pathFile=$outDir/$outPrefix.path_local

#------------------------------------------------------------------------------
if [ "$rootseq" = "" ]; then
  ../bin/epievo_sim $paramFile -t $treeFile \
    -v -n $sites -o $statesFile -p $globalPathFile
else
  python ../pyscripts/write_rootseq_to_file.py $rootseq $rootFile

  ../bin/epievo_sim $paramFile -t $treeFile -r $rootFile \
    -v -n $sites -o $statesFile -p $globalPathFile
fi

../bin/global_jumps_to_paths $treeFile $statesFile $globalPathFile $pathFile
