evosim test.param -o test.paths &> test.out
awk 'NR>6{print}' test.out.tmp
Rscript plot.R test.paths test_evo_path.pdf test.out test_freq_trace.pdf
