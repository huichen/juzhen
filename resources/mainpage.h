/** @mainpage Matrix Library for C++ (mlcpp) 

@author Hui Chen "<usa.chen@gmail.com>" 

@section download Download 

Download the <a href="https://github.com/huichen/mlcpp/tarball/master">tarball</a> or pull the source from <a href="https://github.com/huichen/mlcpp">github</a>:

- git clone https://github.com/huichen/mlcpp.git

@section howtostart How to start

Read the \ref how_to_start page.

@section benchmarks Benchmarks

Speed is the major concern when we design the library. The code has been extensively optimized to make sure it compares favorably to other C++ matrix libraries in terms of speed. 

We have benchmarked our library and then compared the results with <a href="http://eigen.tuxfamily.org">Eigen</a> which provides similar functionalities. The results have shown that mlcpp performs at least equally good as Eigen in most scenarios while it runs significantly faster than Eigen in some key operations such as element indexing and matrix-matrix multiplication.

See the \ref benchmark page for details.

@section problems_and_improvements Problems and improvements

Send me an email or report issues through github:

https://github.com/huichen/mlcpp/issues

I don't have time to write the new features that you want but I'm very glad to merge them into my branch.

*/