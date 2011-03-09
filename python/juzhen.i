%module juzhen
%{
#include <juzhen.h>
typedef juzhen::Complex<double> Complex; 
typedef juzhen::Matrix<double> Matrix; 
typedef juzhen::Matrix<juzhen::Complex<double> > CMatrix; 
typedef juzhen::Vector<double> Vector; 
typedef juzhen::Vector<juzhen::Complex<double> > CVector; 

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

using juzhen::RandVector;
using juzhen::RandMatrix;
%}

%include "complex.i"

%include "matrix.i"

%include "vector.i"

%include "extension.i"

%include "util.i"

