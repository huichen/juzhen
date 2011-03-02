%module mlpy
%{
#include <mlcpp.h>
typedef mlcpp::Complex<double> Complex; 
typedef mlcpp::Matrix<double> Matrix; 
typedef mlcpp::Matrix<mlcpp::Complex<double> > CMatrix; 
%}

%include "mlpy_complex.i"

%include "mlpy_matrix.i"
