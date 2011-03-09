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

#ifndef SRC_UTIL_MAP_REDUCE_H_
#define SRC_UTIL_MAP_REDUCE_H_
#include <math.h>

#include <core/matrix.h>
#include <core/vector.h>

#include <algorithm>
#include <vector>
#include <string>

namespace mlcpp {
/**
 * Map a function func to all elements in a matrix.
 */
template<typename T>
Matrix<T> Map(const Matrix<T> &m, T (*func)(T)) {
  Matrix<T> matrix(m);
  size_t endi = matrix.num_row() * matrix.num_col();
  for (size_t i = 0; i < endi; i++)
    matrix(i) = func(matrix(i));
  matrix.set_temporary(true);
  return matrix;
}

/**
 * Map a function func to all elements in a vector.
 */
template<typename T>
Vector<T> Map(const Vector<T> &v, T (*func)(T)) {
  Vector<T> vector(v);
  size_t endi = vector.size();
  for (size_t i = 0; i < endi; i++)
    vector(i) = func(vector(i));
  vector.set_temporary(true);
  return vector;
}
}
#endif  // SRC_UTIL_MAP_REDUCE_H_
