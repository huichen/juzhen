/*
+---------------------------------------------------------------------------+
|  Juzhen: C++ library for linear algebra                                   |
+---------------------------------------------------------------------------+
|                                                                           |
|  Copyright 2011 Hui Chen                                                  |
|                                                                           |
|  Licensed under the Apache License, Version 2.0 (the "License");          |
|  you may not use this file except in compliance with the License.         |
|  You may obtain a copy of the License at                                  |
|                                                                           |
|      http://www.apache.org/licenses/LICENSE-2.0                           |
|                                                                           |
|  Unless required by applicable law or agreed to in writing, software      |
|  distributed under the License is distributed on an "AS IS" BASIS,        |
|  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. |
|  See the License for the specific language governing permissions and      |
|  limitations under the License.                                           |
|                                                                           |
+---------------------------------------------------------------------------+
*/

#ifndef SRC_UTIL_MATRIX_MATH_H_
#define SRC_UTIL_MATRIX_MATH_H_
#include <core/matrix.h>

namespace juzhen {
/**
 * Return the trace of a matrix 
 */
template<typename T>
T Trace(const Matrix<T> &matrix) {
  assert(matrix.num_col() == matrix.num_row());
  size_t endi = matrix.num_col();
  if (endi == 0)
    return 0;
  T sum = 0;
  for (size_t i = 0; i < endi; i++)
    sum += matrix(i, i);
  return sum;
}

/**
 * scalar + matrix 
 */
template<typename DataType>
inline Matrix<DataType> operator+(
    double lhs,
    const Matrix<DataType> &ma) {
  return ma + lhs;
}

template<typename DataType>
inline Matrix<DataType> operator+(
    const CD &lhs,
    const Matrix<DataType> &ma) {
  return ma + lhs;
}

/**
 * scalar - matrix 
 */
template<typename DataType>
inline Matrix<DataType> operator-(
    double lhs,
    const Matrix<DataType> &ma) {
  return (-ma) + lhs;
}

template<typename DataType>
inline Matrix<DataType> operator-(
    const CD &lhs,
    const Matrix<DataType> &ma) {
  return (-ma) + lhs;
}

/**
 * Multiply a double number and a Matrix
 */
template<typename DataType>
inline Matrix<DataType> operator*(
    double lhs,
    const Matrix<DataType> &ma) {
  return ma*lhs;
}

/**
 * Multiply a complex number and a Matrix
 */
template<typename DataType>
inline Matrix<DataType> operator*(
    const CD &lhs,
    const Matrix<DataType> &ma) {
  return ma*lhs;
}

/**
 * Schur (Hadamard) product of two matrices.
 */
template<typename T>
Matrix<T> SchurProduct(
    const Matrix<T> &matrix1,
    const Matrix<T> &matrix2) {
  assert(matrix1.num_col() == matrix2.num_col());
  assert(matrix1.num_row() == matrix2.num_row());
  Matrix<T> matrix_return(matrix1);

  size_t endi = matrix1.size();
  for (size_t i = 0; i < endi; i++)
    matrix_return(i) *= matrix2(i);
  matrix_return.set_temporary(true);
  return matrix_return;
}

/**
 * Kronecker product of two matrices.
 */
template<typename T>
Matrix<T> KroneckerProduct(
    const Matrix<T> &matrix1,
    const Matrix<T> &matrix2) {
  Matrix<T> matrix(matrix1.num_row() * matrix2.num_row(),
                   matrix1.num_col() * matrix2.num_col());

  size_t endr1, endr2, endc1, endc2;
  endr1 = matrix1.num_row();
  endr2 = matrix2.num_row();
  endc1 = matrix1.num_col();
  endc2 = matrix2.num_col();
  for (size_t r1 = 0; r1 < endr1; r1++)
    for (size_t r2 = 0; r2 < endr2; r2++)
      for (size_t c1 = 0; c1 < endc1; c1++)
        for (size_t c2 = 0; c2 < endc2; c2++)
          matrix(r1 * endr2 + r2, c1 * endc2 + c2) =
              matrix1(r1, c1) * matrix2(r2, c2);
  matrix.set_temporary(true);
  return matrix;
}

/**
 * Inner (Frobenius) product of two matrices.
 * 
 * Inner product of A, B = trace(B* * A)
 */
template<typename T>
T InnerProduct(
    const Matrix<T> &matrix1,
    const Matrix<T> &matrix2) {
  assert(matrix1.num_col() == matrix2.num_col());
  assert(matrix1.num_row() == matrix2.num_row());
  return Trace(Adjoint(matrix2)*matrix1);
}

template<typename T>
T InnerProduct(const Matrix<T> &matrix) {
  return InnerProduct(matrix, matrix);
}

template<typename T>
T InnerProduct(const Matrix<Complex<T> > &matrix) {
  return InnerProduct(matrix, matrix).real;
}
}
#endif  // SRC_UTIL_MATRIX_MATH_H_
