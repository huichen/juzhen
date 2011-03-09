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

#ifndef SRC_UTIL_RANDOM_H_
#define SRC_UTIL_RANDOM_H_
#include <stdlib.h>
#include <time.h>

#include <core/matrix.h>
#include <core/vector.h>

namespace mlcpp {
/**
 * Return a random Vector with elements between 0 and scale.
 */
template<typename T, typename T1>
Vector<T> &Randomize(Vector<T> &vector, T1 scale) {  // NOLINT
  unsigned int seed = (unsigned int)time(0);
  size_t size = vector.size();
  for (size_t i = 0; i < size; i++)
    vector(i) = scale * (T)rand_r(&seed) / (T)RAND_MAX;
  return vector;
}

/**
 * Return a random Matrix with elements between 0 and scale.
 */
template<typename T, typename T1>
Matrix<T>& Randomize(Matrix<T> &matrix, T1 scale) {  // NOLINT
  unsigned int seed = (unsigned int)time(0);
  size_t size = matrix.num_col() * matrix.num_row();
  for (size_t i = 0; i < size; i++)
    matrix(i) = scale * (T)rand_r(&seed) / (T)RAND_MAX;
  return matrix;
}
}
#endif  // SRC_UTIL_RANDOM_H_
