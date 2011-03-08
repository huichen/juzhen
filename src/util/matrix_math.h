/*
+---------------------------------------------------------------------------+
|  Matrix Library for C++ (mlcpp)                                           |
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

namespace mlcpp {
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
  for (size_t i = 1; i < endi; i++)
    sum += matrix(i, i);
  return sum;
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

  size_t endi = matrix1.num_col() * matrix1.num_row();
  for (size_t i = 0; i < endi; i++)
    matrix_return(i) *= matrix2(i);
  matrix_return.set_temporary(true);
  return matrix_return;
}

/**
 * Inner product of two matrices.
 * 
 * Inner product of A, B = trace(B* * A)
 */
template<typename T>
T InnerProduct(
    const Matrix<T> &matrix1,
    const Matrix<T> &matrix2) {
  assert(matrix1.num_col() == matrix2.num_col());
  assert(matrix1.num_row() == matrix2.num_row());
  return trace(Adjoint(matrix2)*matrix1);
}

template<typename T>
T InnerProduct(const Matrix<T> &matrix) {
  return InnerProduce(matrix, matrix);
}

template<typename T>
T InnerProduct(const Matrix<Complex<T> > &matrix) {
  return InnerProduce(matrix, matrix).real;
}
}
#endif  // SRC_UTIL_MATRIX_MATH_H_
