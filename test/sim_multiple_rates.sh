number=100
sites=1000
time=10.0
input="input/rates/multiple_rates.txt"
output=multiple_rates

print_usage() {
  printf "Usage: $(basename $0) [-n number] [-s sites] [-t time]
         [-a T00] [-b T11] [-i input rates] [-o output] [-h help] \n"
}

while getopts 'n:s:t:i:o:Rh' flag; do
  case "${flag}" in
    n) number="${OPTARG}" ;;
    s) sites="${OPTARG}" ;;
    t) time="${OPTARG}" ;;
    i) input="${OPTARG}" ;;
    o) output="${OPTARG}" ;;
    R) sampleRoot=1 ;;
    h) print_usage
       exit 0 ;;
    *) print_usage
       exit 1 ;;
  esac
done

inputDir=input
outputDir=output
paramDir=$inputDir/rates
mkdir -p $paramDir

while IFS='' read -r line || [[ -n "$line" ]]; do
  #############################################################################
  # WRITING INPUT FILES
  #############################################################################
  B=`echo $line | cut -d' ' -f1`
  D=`echo $line | cut -d' ' -f2`
  E=`echo $line | cut -d' ' -f3`
  C=`echo $line | cut -d' ' -f4`
  M=`echo $line | cut -d' ' -f5`
  init0=`echo $line | cut -d' ' -f6`
  init1=`echo $line | cut -d' ' -f7`
  S=$(python <( echo 'import sys;
B=float(sys.argv[1]);
C=float(sys.argv[2]);
M=float(sys.argv[3]);
E=float(sys.argv[4]);
D=float(sys.argv[5]);
print((B * C * C * M) / (E * E * D))' ) $B $C $M $E $D)

  ### Generate parameter file
  prefix=${B}_${D}_${E}_${C}_${M}_${init0}_${init1}
  paramFile=$paramDir/$prefix.rates.tmp
  echo -e "000\t"$B"\n001\t"$E"\n010\t"$D"\n011\t"$C \
    "\n100\t"$E"\n101\t"$M"\n110\t"$C"\n111\t"$S > $paramFile

  ### Generate tree file
  treeFile=$inputDir/tree_T${time}.nwk
  echo -e "(C:"$time",D:0.0)E:0.0;" > $treeFile

  ### Generate root sequence
  rootFile=$inputDir/rootseq/${init0}_${init1}.rootseq
  if [ ! -f $rootFile ];then
    ./sample_root.sh -s $sites -a $init0 -b $init1 -o $rootFile
  fi
  #############################################################################
  # SIMULATING
  #############################################################################
  echo "Simulating: " $prefix

  pathDir=$outputDir/$output/paths
  statDir=$outputDir/$output/stats
  mkdir -p $pathDir
  mkdir -p $statDir

  stateFile=$pathDir/$prefix.outmeth
  pathFile=$pathDir/$prefix.global_jumps
  statFile=$statDir/$prefix.stats
  echo -e "A_00\tA_01\tA_10\tA_11\tN_B\tN_D\tN_E\tN_C\tN_M\tN_S" \
    "\tA_B\tA_D\tA_E\tA_C\tA_M\tA_S" > $statFile

  for (( i=0; i<$number; i++ )); do
    epievo_sim $paramFile -r $rootFile -T $time -o $stateFile -n $sites -p $pathFile
    evo_features.py -i $stateFile |tail -n 1 >> $statFile
  done

done < $input

###############################################################################
# CLEAN UP
###############################################################################
rm $paramDir/*.tmp
