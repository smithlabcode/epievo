num=1
sites=5
prefix=""
rootseq=""
dpar=""
dsim=$PWD

print_usage() {
  printf "Usage: $(basename $0) [-n number] [-s sites] [-r root sequence]
          [-f file prefix] [-p parameter directory] [-o output directory]\n"
}

while getopts 'n:s:r:f:p:o:h' flag; do
  case "${flag}" in
    n) num="${OPTARG}" ;;
    s) sites="${OPTARG}" ;;
    r) rootseq="${OPTARG}" ;;
    f) prefix="${OPTARG}" ;;
    p) dpar="${OPTARG}" ;;
    o) dsim="${OPTARG}" ;;
    h) print_usage
       exit 0 ;;
    *) print_usage
       exit 1 ;;
  esac
done

outprefix=${prefix}_N${sites}

param=$dpar/$prefix.param
tree=$dpar/$prefix.nwk
rootfile=$dpar/$rootseq.rootseq

outmeth=$dsim/$outprefix.outmeth
jumps=$dsim/$outprefix.global_jumps
paths=$dsim/$outprefix.path_local
stats=$dsim/$outprefix.stats

#------------------------------------------------------------------------------
if [ "$rootseq" = "" ]; then
  stationary_stats_est $param -t $tree \
    -n $sites -m $num \
    -o $outmeth -p $jumps -S $stats
else
  write_rootseq_to_file.py $rootseq $rootfile

  stationary_stats_est $param -t $tree -r $rootfile \
    -n $sites -m $num \
    -o $outmeth -p $jumps -S $stats
fi

global_jumps_to_paths $tree $outmeth $jumps $paths
