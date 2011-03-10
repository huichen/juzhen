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

#ifndef SRC_UTIL_VECTOR_FITTING_H_
#define SRC_UTIL_VECTOR_FITTING_H_
#include <math.h>

#include <core/vector.h>

namespace juzhen {
/**
 * Least square method.
 * Find the factors k1 and k2 to make sum(|k1*a + k2 - b|^2) minimum.
 */
template<typename T>
Vector<T> LinearLeastSquares(const Vector<T> &a, const Vector<T> &b) {
  assert(a.size() == b.size());
  int n = a.size();
  T mean_a = Mean(a);
  T mean_b = Mean(b);
  Vector<T> vector(2);
  vector.Clear();

  T temp = a * a - n * mean_a * mean_a;
  if (temp != 0) {
    vector(0) = (a * b - n * mean_a * mean_b) / temp;
    vector(1) = mean_b - mean_a * vector(0);
  }

  vector.set_temporary(true);
  return vector;
}
}
#endif  // SRC_UTIL_VECTOR_FITTING_H_
