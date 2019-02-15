rep=$1

dsrc=/Users/lizji/codes/EvoSim/src
prefix=test
dout=output/forward

cd /Users/lizji/data/evo/mle

> $dout/st0.stats
> $dout/st1.stats
> $dout/bl0.stats
> $dout/bl1.stats
> $dout/rt0.stats
> $dout/rt1.stats

for (( i=0; i<$rep; i++ )); do
  $dsrc/prog/epievo_sim input/${prefix}.param -t input/tree.nwk -o input/${prefix}.outmeth -n 10000 -p input/${prefix}.global_jumps
  $dsrc/prog/global_jumps_to_paths input/tree.nwk input/${prefix}.outmeth input/${prefix}.global_jumps input/${prefix}.path_local
  $dsrc/prog/epievo_est_complete -S -p input/${prefix}.param -t input/tree.nwk -o $dout/${prefix}.param.update input/${prefix}.path_local
 
  paramFile=output/${prefix}.param.update
  text=$(<$paramFile)
  st0=`echo $text | cut -d' ' -f2`
  st1=`echo $text | cut -d' ' -f3`
  bl0=`echo $text | cut -d' ' -f5`
  bl1=`echo $text | cut -d' ' -f6`
  rt0=`echo $text | cut -d' ' -f8`
  rt1=`echo $text | cut -d' ' -f9`

  echo $st0 >> $dout/st0.stats
  echo $st1 >> $dout/st1.stats
  echo $bl0 >> $dout/bl0.stats
  echo $bl1 >> $dout/bl1.stats
  echo $rt0 >> $dout/rt0.stats
  echo $rt1 >> $dout/rt1.stats

done

/Users/lizji/data/evo/scripts/mle/mle_forward.R
