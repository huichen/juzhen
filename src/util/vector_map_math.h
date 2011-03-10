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
#ifndef SRC_UTIL_VECTOR_MAP_MATH_H_
#define SRC_UTIL_VECTOR_MAP_MATH_H_
#include <math.h>

#include <core/vector.h>
#include <util/map_reduce.h>

namespace juzhen {

/// Trigonometric functions
template<typename T>
Vector<T> Cos(const Vector<T> &vector) {
  return Map(vector, &cos);
}

template<typename T>
Vector<T> Sin(const Vector<T> &vector) {
  return Map(vector, &sin);
}

template<typename T>
Vector<T> Tan(const Vector<T> &vector) {
  return Map(vector, &tan);
}

template<typename T>
Vector<T> Acos(const Vector<T> &vector) {
  return Map(vector, &acos);
}

template<typename T>
Vector<T> Asin(const Vector<T> &vector) {
  return Map(vector, &asin);
}

template<typename T>
Vector<T> Atan(const Vector<T> &vector) {
  return Map(vector, &atan);
}

/// Hyperbolic functions
template<typename T>
Vector<T> Cosh(const Vector<T> &vector) {
  return Map(vector, &cosh);
}

template<typename T>
Vector<T> Sinh(const Vector<T> &vector) {
  return Map(vector, &sinh);
}

template<typename T>
Vector<T> Tanh(const Vector<T> &vector) {
  return Map(vector, &tanh);
}

/// Exponential and logarithmic functions

template<typename T>
Vector<T> Exp(const Vector<T> &vector) {
  return Map(vector, &exp);
}

template<typename T>
Vector<T> Log(const Vector<T> &vector) {
  return Map(vector, &log);
}

template<typename T>
Vector<T> Log10(const Vector<T> &vector) {
  return Map(vector, &log10);
}

/// Power functions

template<typename T>
Vector<T> Sqrt(const Vector<T> &vector) {
  return Map(vector, &sqrt);
}

/// Rounding, absolute value and remainder functions
template<typename T>
Vector<T> Ceil(const Vector<T> &vector) {
  return Map(vector, &ceil);
}

template<typename T>
Vector<T> Fabs(const Vector<T> &vector) {
  return Map(vector, &fabs);
}

template<typename T>
Vector<T> Floor(const Vector<T> &vector) {
  return Map(vector, &floor);
}
}
#endif  // SRC_UTIL_VECTOR_MAP_MATH_H_
