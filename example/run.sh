prefix=test
evosim ${prefix}.param -v -p ${prefix}_path.txt -o ${prefix}_meth.txt
Rscript plot.R ${prefix}_path.txt ${prefix}_path.pdf 

