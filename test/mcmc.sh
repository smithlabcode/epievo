sites=5
num=10
batch=10
evotime=-1
seed=-1

fixRoot=0
fullPath=0

paramFile=test.param
#paramFile=input/paired.rates
#paramFile=input/indep.rates
treeFile=tree.nwk
outPrefix=test
outDir=output/mcmc

print_usage() {
  printf "Usage: $(basename $0)
          [-n sites ($sites)]
          [-i number of MCMC batches ($num)]
          [-B MCMC batch ($batch)]
          [-T evolutionary time]
          [-s seed ($seed)]
          [-R fix root ($fixRoot)]
          [-F full-path sampling ($fullPath)]
          [-p param file ($paramFile)]
          [-t tree file ($treeFile)]
          [-f output file prefix ($outPrefix)]
          [-o output directory ($outDir)]\n"
}

while getopts 'n:i:B:T:s:RFp:t:f:o:h' flag; do
  case "${flag}" in
    n) sites="${OPTARG}" ;;
    i) num="${OPTARG}" ;;
    B) batch="${OPTARG}" ;;
    T) evotime="${OPTARG}" ;;
    s) seed="${OPTARG}" ;;
    R) fixRoot=1 ;;
    F) fullPath=1 ;;
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

statsFile=$outDir/$outPrefix
figureFile=$outDir/$outPrefix.pdf

#------------------------------------------------------------------------------
echo "1. Sampling"
CMD="../bin/mcmc_test -i $num -n $sites -B $batch
-v -o $statsFile $paramFile"
if [ "$evotime" != -1 ];then
  CMD="$CMD -T $evotime"
else
  CMD="$CMD -t $treeFile"
fi

if [ "$fullPath" -eq 1 ];then
  CMD="$CMD -F"
fi

if [ "$fixRoot" -eq 1 ];then
  CMD="$CMD -R"
fi

if [ "$seed" -gt 0 ];then
  CMD="$CMD -s $seed"
fi

  echo $CMD
eval $CMD

#------------------------------------------------------------------------------
echo "2. Plotting"
Rscript ../rscripts/plot_sufficient_stats.R -p $outPrefix -i $outDir \
  -o $figureFile
