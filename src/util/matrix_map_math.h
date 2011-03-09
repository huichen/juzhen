/*
+---------------------------------------------------------------------------+
|  Juzhen: C++ library for linear algebra                                           |
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
#ifndef SRC_UTIL_MATRIX_MAP_MATH_H_
#define SRC_UTIL_MATRIX_MAP_MATH_H_
#include <math.h>

#include <core/matrix.h>
#include <util/map_reduce.h>

namespace juzhen {

/// Trigonometric functions
template<typename T>
Matrix<T> Cos(const Matrix<T> &matrix) {
  return Map(matrix, &cos);
}

template<typename T>
Matrix<T> Sin(const Matrix<T> &matrix) {
  return Map(matrix, &sin);
}

template<typename T>
Matrix<T> Tan(const Matrix<T> &matrix) {
  return Map(matrix, &tan);
}

template<typename T>
Matrix<T> Acos(const Matrix<T> &matrix) {
  return Map(matrix, &acos);
}

template<typename T>
Matrix<T> Asin(const Matrix<T> &matrix) {
  return Map(matrix, &asin);
}

template<typename T>
Matrix<T> Atan(const Matrix<T> &matrix) {
  return Map(matrix, &atan);
}

/// Hyperbolic functions
template<typename T>
Matrix<T> Cosh(const Matrix<T> &matrix) {
  return Map(matrix, &cosh);
}

template<typename T>
Matrix<T> Sinh(const Matrix<T> &matrix) {
  return Map(matrix, &sinh);
}

template<typename T>
Matrix<T> Tanh(const Matrix<T> &matrix) {
  return Map(matrix, &tanh);
}

/// Exponential and logarithmic functions

template<typename T>
Matrix<T> Exp(const Matrix<T> &matrix) {
  return Map(matrix, &exp);
}

template<typename T>
Matrix<T> Log(const Matrix<T> &matrix) {
  return Map(matrix, &log);
}

template<typename T>
Matrix<T> Log10(const Matrix<T> &matrix) {
  return Map(matrix, &log10);
}

/// Power functions

template<typename T>
Matrix<T> Sqrt(const Matrix<T> &matrix) {
  return Map(matrix, &sqrt);
}

/// Rounding, absolute value and remainder functions
template<typename T>
Matrix<T> Ceil(const Matrix<T> &matrix) {
  return Map(matrix, &ceil);
}

template<typename T>
Matrix<T> Fabs(const Matrix<T> &matrix) {
  return Map(matrix, &fabs);
}

template<typename T>
Matrix<T> Floor(const Matrix<T> &matrix) {
  return Map(matrix, &floor);
}
}
#endif  // SRC_UTIL_MATRIX_MAP_MATH_H_
