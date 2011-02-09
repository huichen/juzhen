/** 
\page how_to_start How to start

@section prerequisite Prerequisite

You need a package that provides CBLAS functions. I recommend optimized BLAS libraries such as GotoBLAS2 and Intel's Math Kernel Library. 

@section usage How to use mlcpp 

-# Set MLCPP_BLASLIB envron to mkl if you want to use MKL or GotoBlas2<br> <br> 
   export MLCPP_BLASLIB=mkl<br><br> 
   or <br><br> 
   export MLCPP_BLASLIB=gotoblas2<br><br> 
-# Run unit tests and make sure it passes everything:<br><br> 
   cd mlcpp_source_dir <br> <br> 
   make <br><br> 
-# Include mlcpp.hpp in your source file and add '-Imlcpp_source_dir' to your compiler flags. <br><br> 
   Use '-DUSE_MKL' if you want to link against MKL.<br><br> 

@section example Examples
  
The package includes a very simple example (example.cpp). There are also plenty sample usages of mlcpp in unittest.cpp file. 

Mlcpp's functions are similar to those of other C++ matrix libraries such as uBlas and Eigen. Read the <a href="annotated.html">API page</a> for reference. 

@section notes Notes

-# Following functions only compile with MKL: <b>matrix::eigen</b>, <b>matrix::leigen</b>, <b>matrix::reigen</b> <br><br> 
-# matrix::getDataPtr() provides column-major raw pointer for low-level BLAS/LAPACK functions.<br><br> 
-# Avoid declaring a matrix or a vector using copy constructor. It will fail. 

*
*/
