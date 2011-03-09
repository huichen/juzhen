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
 * Find the factor k to make |k * a + b| minimum.
 */
template<typename T>
T LeastSquaresMethod(const Vector<T> &a, const Vector<T> &b) {
  return -(a * b) / (a * a);
}
}
#endif  // SRC_UTIL_VECTOR_FITTING_H_
