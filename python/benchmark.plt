set term png 
set out "benchmark_mlcpp.png"
set title "Matrix Product (10 times)" 
set xlabel "Matrix Rank" 
set ylabel "Seconds"
plot 'benchmark.txt' using 1:2 title 'mlcpp' w l, 'benchmark.txt' using 1:3 title 'Numpy' w l 

