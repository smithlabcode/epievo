prefix=test
evosim ${prefix}.param -v -o ${prefix}_path.txt 1> ${prefix}_freq.txt 
Rscript plot.R ${prefix}_path.txt ${prefix}_path.pdf ${prefix}_freq.txt  ${prefix}_freq.pdf 

