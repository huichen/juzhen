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

#ifndef SRC_UTIL_MATRIX_STATISTICS_H_
#define SRC_UTIL_MATRIX_STATISTICS_H_
#include <core/matrix.h>

namespace juzhen {
/**
 * Return the maximum element of a matrix 
 */
template<typename T>
T Max(const Matrix<T> &matrix) {
  size_t endi = matrix.size();
  if (endi == 0)
    return 0;
  T max_value = matrix(0);
  for (size_t i = 1; i < endi; i++)
    if (max_value < matrix(i)) max_value = matrix(i);
  return max_value;
}

/**
 * Return the minimum element of a matrix 
 */
template<typename T>
T Min(const Matrix<T> &matrix) {
  size_t endi = matrix.size();
  if (endi == 0)
    return 0;
  T min_value = matrix(0);
  for (size_t i = 1; i < endi; i++)
    if (min_value > matrix(i)) min_value = matrix(i);
  return min_value;
}

/**
 * Return the sum of all elements in a matrix 
 */
template<typename T>
T Sum(const Matrix<T> &matrix) {
  size_t endi = matrix.size();
  if (endi == 0)
    return 0;
  T sum = 0;
  for (size_t i = 0; i < endi; i++)
    sum += matrix(i);
  return sum;
}

/**
 * Return the product of all elements in a matrix 
 */
template<typename T>
T Prod(const Matrix<T> &matrix) {
  size_t endi = matrix.size();
  if (endi == 0)
    return 0;
  T prod = 1;
  for (size_t i = 0; i < endi; i++)
    prod *= matrix(i);
  return prod;
}

/**
 * Return the average of all elements in a matrix 
 */
template<typename T>
T Mean(const Matrix<T> &matrix) {
  size_t endi = matrix.size();
  if (endi == 0)
    return 0;
  T sum = 0;
  for (size_t i = 0; i < endi; i++)
    sum += matrix(i);
  return sum/endi;
}

/**
 * Return the norm-square of a Matrix.
 */
template<typename T>
double NormSquare(const Matrix<T> &m) {
  double r = 0;
  size_t size = m.size();
  for (size_t i = 0; i < size; i++)
    r += abs2(m(i));
  return r;
}

/**
 * Return the norm of a Matrix.
 */
template<typename T>
double Norm(const Matrix<T> &m) {
  return sqrt(NormSquare(m));
}

/**
 * Return the 1-norm of a Matrix.
 */
template<typename T>
double NormOne(const Matrix<T> &m) {
  double r = 0;
  size_t size = m.size();
  for (size_t i = 0; i < size; i++)
    r += abs(m(i));
  return r;
}

/**
 * Return the infinity-norm of a Matrix.
 */
template<typename T>
double NormInfinity(const Matrix<T> &m) {
  size_t size = m.size();
  if (size == 0)
    return 0;
  double r = abs(m(0));
  for (size_t i = 1; i < size; i++)
    if (r < abs(m(i)))
      r = abs(m(i));
  return r;
}
}
#endif  // SRC_UTIL_MATRIX_STATISTICS_H_
