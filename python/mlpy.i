%module mlpy
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
%}

%include "mlpy_complex.i"

%include "mlpy_matrix.i"

%include "mlpy_vector.i"
