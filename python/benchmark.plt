set term png 
set out "benchmark_juzhen.png"
set title "Matrix Product (10 times)" 
set xlabel "Matrix Size" 
set ylabel "Seconds"
plot 'benchmark.txt' using 1:2 title 'Juzhen' w l, 'benchmark.txt' using 1:3 title 'Numpy' w l 

