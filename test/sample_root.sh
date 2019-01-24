sites=1000
init0=0.95
init1=0.9
verbose=0

dout=input/rootseq
mkdir -p $dout
output=$dout/${init0}_${init1}.rootseq

print_usage() {
  printf "Usage: $(basename $0) [-s sites] [-a T00] [-b T11] [-o output]
          [-v verbose] \n"
}

while getopts 's:a:b:o:v' flag; do
  case "${flag}" in
    s) sites="${OPTARG}" ;;
    a) init0="${OPTARG}" ;;
    b) init1="${OPTARG}" ;;
    o) output="${OPTARG}" ;;
    v) verbose=1 ;;
    h) print_usage
       exit 0 ;;
    *) print_usage
       exit 1 ;;
  esac
done

#------------------------------------------------------------------------------
#echo "------------------------------------------------------------------------"
#echo Number of sites: $sites
echo T00: $init0
echo T11: $init1
#echo Output: $output
########### STEP 1 ############################################################
echo "------------------------------------------------------------------------"
echo 1. Write parameter file...
fparam=$output.param.tmp
st0=$init0
st1=$init1
bl0=-1
bl1=-1

echo -e "stationary\t"$st0"\t"$st1"\nbaseline\t"$bl0"\t"$bl1"\ninit\t"$init0"\t"$init1 > $fparam

########### STEP 2 ############################################################
echo Simulating...
simCMD="epievo_sim $fparam -T 0.1 -o $output -n $sites -p ${output}.global_jumps.tmp"
if [ "$verbose" -eq 1 ]; then
  simCMD="${simCMD} -v"
fi
eval $simCMD

########### STEP 3 ############################################################
echo  Clean-up...

rm $output.*.tmp
