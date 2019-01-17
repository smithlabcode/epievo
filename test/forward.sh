num=1000
sites=5
prefix=""
dpar=""
dsim=$PWD
dfw=$PWD
fixRoot=0

print_usage() {
  printf "Usage: $(basename $0) [-n number] [-s sites] [-f file prefix]
          [-f file prefix] [-p parameter directory]
          [-i input path directory] [-o output directory] [-R fix root]\n"
}

while getopts 'n:s:f:p:i:o:Rh' flag; do
  case "${flag}" in
    n) num="${OPTARG}" ;;
    s) sites="${OPTARG}" ;;
    f) prefix="${OPTARG}" ;;
    p) dpar="${OPTARG}" ;;
    i) dsim="${OPTARG}" ;;
    o) dfw="${OPTARG}" ;;
    R) fixRoot=1 ;;
    h) print_usage
       exit 0 ;;
    *) print_usage
       exit 1 ;;
  esac
done

outprefix=${prefix}_N${sites}

param=$dpar/$prefix.param
tree=$dpar/$prefix.nwk
inpaths=$dsim/$outprefix.path_local
outpaths=$dfw/$outprefix.path_local
stats=$dfw/$outprefix.stats

#------------------------------------------------------------------------------
if [ "$fixRoot" -eq 1 ]; then
  forward_prop_test -r $num -S $stats -o $outpaths \
    -R $param $tree $inpaths
else
  forward_prop_test -r $num -S $stats -o $outpaths \
    $param $tree $inpaths
fi
