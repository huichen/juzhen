/** @mainpage Matrix Library for C++ (mlcpp) 
*
* @author Hui Chen "<usa dot chen at gmail>" 
*

* @section prerequisite Prerequisite
  You need a package that provides CBLAS functions, such as MKL. 

* @section download Download 
  
  Pull the source code of mlcpp from github:

  https://github.com/huichen/mlcpp

  Direct link to tarball is

  https://github.com/huichen/mlcpp/tarball/master 

* @section how_to_use How to use
  -# Set MLCPP_BLASLIB envron to mkl if you want to use MKL <br> <br> 
     export MLCPP_BLASLIB=mkl<br><br> 
  -# Run unit tests and make sure it passes everything:<br><br> 
     cd mlcpp_source_dir <br> <br> 
     make <br><br> 
  -# Include mlcpp.hpp in your source file and add '-Imlcpp_source_dir' to your compiler flags. <br><br> 
     Use '-DUSE_MKL' if you build with MKL <br><br> 

* @section example Examples
  
  There are plenty samples in unittest.cpp file. Reading comments in the code helps too.

* @section notes Notes

  -# Following functions only compile with MKL <br><br> 
     matrix::eigen, matrix::leigen, matrix::reigen <br><br> 
  -# matrix::getDataPtr() provides column-major raw pointer for low-level BLAS/LAPACK functions.<br><br> 

*  @section problems_and_improvements Problems and improvements

  Send me an email or report issues through github:

  https://github.com/huichen/mlcpp/issues

  I don't have time to write the new features that you want but I'm very glad to merge them into my branch.

*
*/
