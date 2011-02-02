/** @mainpage package templates
*
* @author Hui Chen "<usa dot chen at gmail>" 
*

* @section prerequisite Prerequisite
  You need a package that provides CBLAS functions, such as MKL. 

* @section how_to_use How to use
  -# Set MLCPP_BLASLIB envron to mkl if you want to use MKL <br> <br> 
     export MLCPP_BLASLIB=mkl<br><br> 
  -# Run unit tests first and make sure it passes everything:<br><br> 
     cd mlcpp_source_dir<br><br> 
     make <br><br> 
  -# #include <mlcpp.hpp> in your source file <br><br> 
     and add '-Imlcpp_source_dir' to your compile flags. <br><br> 
     Use '-DUSE_MKL' if you use MKL <br><br> 

* @section notes Notes

  -# Following functions only compiles with MKL <br><br> 
     matrix::eigen, matrix::leigen, matrix::reigen <br><br> 
  -# matrix::getDataPtr() provides column-major raw array pointer of a matrix to use in low-level BLAS/LAPACK functions  <br><br> 

*  @section problems_and_improvements Problems and improvements

  Please report issues through github:

  https://github.com/huichen/mlcpp/issues

  I don't have time to write all new features you want but I'm very glad to merge them into my branch.

*
*/
