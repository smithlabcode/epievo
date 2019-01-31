num=1000
sites=5
proposal=2
estParam=0
sampleTree=0
fixRoot=0
burning=1000
batch=1
prefix=""
initprefix=$prefix
dpar=""
dsim=$PWD
dmcmc=$PWD

print_usage() {
  printf "Usage: $(basename $0) [-n number] [-s sites] [-P proposal] \
          [-E estimate parameters] [-T sample full tree] [-R fix root]
          [-L burning length] [-B batch size]
          [-f file prefix] [-j init parameter prefix] [-p parameter directory]
          [-i input path directory] [-o output directory]\n"
}

while getopts 'n:s:f:j:P:ETRL:B:p:i:o:h' flag; do
  case "${flag}" in
    n) num="${OPTARG}" ;;
    s) sites="${OPTARG}" ;;
    f) prefix="${OPTARG}"
       initprefix=$prefix ;;
    j) initprefix="${OPTARG}" ;;
    P) proposal="${OPTARG}" ;;
    E) estParam=1 ;;
    T) sampleTree=1 ;;
    R) fixRoot=1 ;;
    L) burning="${OPTARG}" ;;
    B) batch="${OPTARG}" ;;
    p) dpar="${OPTARG}" ;;
    i) dsim="${OPTARG}" ;;
    o) dmcmc="${OPTARG}" ;;
    h) print_usage
       exit 0 ;;
    *) print_usage
       exit 1 ;;
  esac
done


outprefix=${prefix}_N${sites}
param=$dpar/$initprefix.param
tree=$dpar/$initprefix.nwk
inpaths=$dsim/$outprefix.path_local
outpaths=$dmcmc/$outprefix.update.path_local
stats=$dmcmc/$outprefix.stats
trace=$dmcmc/$outprefix.trace

#------------------------------------------------------------------------------
mcmcCMD="metropolis_unif_test -P $proposal -B $batch -r $num -L $burning
-S $stats -t $trace -o $outpaths $param $tree $inpaths"

if [ "$estParam" -eq 1 ]; then
  mcmcCMD="${mcmcCMD} -m"
fi
if [ "$fixRoot" -eq 1 ]; then
  mcmcCMD="${mcmcCMD} -R"
fi
if [ "$sampleTree" -eq 0 ]; then
  mcmcCMD="${mcmcCMD} -b 1"
fi

eval $mcmcCMD
