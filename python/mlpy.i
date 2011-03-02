%module mlpy
%{
#include <mlcpp.h>
typedef mlcpp::Complex<double> Complex; 
typedef mlcpp::Matrix<double> Matrix; 
typedef mlcpp::Matrix<mlcpp::Complex<double> > CMatrix; 
typedef mlcpp::Vector<double> Vector; 
typedef mlcpp::Vector<mlcpp::Complex<double> > CVector; 
%}

%include "mlpy_complex.i"

%include "mlpy_matrix.i"

%include "mlpy_vector.i"
