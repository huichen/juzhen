/** 
\page how_to_start How to start

@section usage Usage

  -# Set MLCPP_BLASLIB envron to mkl if you want to use MKL <br> <br> 
     export MLCPP_BLASLIB=mkl<br><br> 
  -# Run unit tests and make sure it passes everything:<br><br> 
     cd mlcpp_source_dir <br> <br> 
     make <br><br> 
  -# Include mlcpp.hpp in your source file and add '-Imlcpp_source_dir' to your compiler flags. <br><br> 
     Use '-DUSE_MKL' if you build with MKL <br><br> 

@section example Examples
  
  There are plenty samples in unittest.cpp file. Reading comments in the code helps too.

@section notes Notes

  -# Following functions only compile with MKL <br><br> 
     matrix::eigen, matrix::leigen, matrix::reigen <br><br> 
  -# matrix::getDataPtr() provides column-major raw pointer for low-level BLAS/LAPACK functions.<br><br> 

*
*/
