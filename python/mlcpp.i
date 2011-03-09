%module mlcpp
%{
#include <mlcpp.h>
typedef mlcpp::Complex<double> Complex; 
typedef mlcpp::Matrix<double> Matrix; 
typedef mlcpp::Matrix<mlcpp::Complex<double> > CMatrix; 
typedef mlcpp::Vector<double> Vector; 
typedef mlcpp::Vector<mlcpp::Complex<double> > CVector; 

Matrix Identity(int i) {
  Matrix m = Matrix(i, i);
  m.Clear();
  for (int j = 0; j < i; j++)
    m(j, j) = 1;
  return m;
}

CMatrix CIdentity(int i) {
  CMatrix m = CMatrix(i, i);
  m.Clear();
  for (int j = 0; j < i; j++)
    m(j, j) = 1;
  return m;
}

using mlcpp::RandVector;
using mlcpp::RandMatrix;
%}

%include "mlcpp_complex.i"

%include "mlcpp_matrix.i"

%include "mlcpp_vector.i"

%include "mlcpp_extension.i"

%include "mlcpp_util.i"

