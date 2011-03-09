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

#ifndef SRC_UTIL_RANDOM_H_
#define SRC_UTIL_RANDOM_H_
#include <stdlib.h>
#include <time.h>

#include <core/matrix.h>
#include <core/vector.h>

namespace juzhen {

static unsigned int sRandomSeed = 9527;

void UpdateRandomSeed() {
  sRandomSeed = (unsigned int)time(0);
}

/**
 * Return a random Vector with elements between 0 and scale.
 */
template<typename T, typename T1>
Vector<T> &Randomize(Vector<T> &vector, T1 scale) {  // NOLINT
  size_t size = vector.size();
  for (size_t i = 0; i < size; i++)
    vector(i) = scale * (T)rand_r(&sRandomSeed) / (T)RAND_MAX;
  return vector;
}

template<typename T>
Vector<T> &Randomize(Vector<T> &vector) {  // NOLINT
  return Randomize(vector, 1);
}

/**
 * Return a random Matrix with elements between 0 and scale.
 */
template<typename T, typename T1>
Matrix<T>& Randomize(Matrix<T> &matrix, T1 scale) {  // NOLINT
  size_t size = matrix.size();
  for (size_t i = 0; i < size; i++)
    matrix(i) = scale * (T)rand_r(&sRandomSeed) / (T)RAND_MAX;
  return matrix;
}

template<typename T>
Matrix<T>& Randomize(Matrix<T> &matrix) {  // NOLINT
  return Randomize(matrix, 1);
}

/**
 * Return a random Vector with elements between 0 and scale.
 */
template<typename T>
Vector<T> RandVector(size_t size, T scale) {  // NOLINT
  Vector<T> vector(size);
  Randomize(vector, scale);
  vector.set_temporary(true);
  return vector;
}

template<typename T>
Vector<T> RandVector(size_t size) {  // NOLINT
  return RandVector<T>(size, 1);
}

/**
 * Return a random Matrix with elements between 0 and scale.
 */
template<typename T>
Matrix<T> RandMatrix(size_t row, size_t col, T scale) {  // NOLINT
  Matrix<T> matrix(row, col);
  Randomize(matrix, scale);
  matrix.set_temporary(true);
  return matrix;
}

template<typename T>
Matrix<T> RandMatrix(size_t row, size_t col) {  // NOLINT
  return RandMatrix<T>(row, col, 1);
}
}
#endif  // SRC_UTIL_RANDOM_H_
