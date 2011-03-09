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

#ifndef SRC_UTIL_VECTOR_BASIC_H_
#define SRC_UTIL_VECTOR_BASIC_H_
#include <math.h>

#include <core/vector.h>

#include <algorithm>
#include <vector>
#include <string>

namespace mlcpp {
/**
 * Return the transpose of a Vector.
 */
template<typename T>
Vector<T> Transpose(const Vector<T> &v) {
  return v;
}

/**
 * Return the hermitian of a Vector.
 */
template<typename T>
Vector<T> Adjoint(const Vector<T> &v) {
  return v.Adjoint();
}

/**
 * Return the conjugate of a Vector.
 */
template<typename T>
Vector<T> Conjugate(const Vector<T> &v) {
  return v.Conjugate();
}

/**
 * Return the real part of a Vector.
 */
Vector<float> Real(const Vector<float> &v) {
  return Real((Matrix<float>)v);
}

/**
 * Return the imaginary part of a Vector.
 */
Vector<float> Imag(const Vector<float> &v) {
  return Imag((Matrix<float>)v);
}

/**
 * Return the real part of a Vector.
 */
Vector<double> Real(const Vector<double> &v) {
  return Real((Matrix<double>)v);
}

/**
 * Return the imaginary part of a Vector.
 */
Vector<double> Imag(const Vector<double> &v) {
  return Imag((Matrix<double>)v);
}

/**
 * Return the real part of a Vector.
 */
Vector<float> Real(const Vector<CS> &v) {
  return Real((Matrix<CS>)v);
}

/**
 * Return the imaginary part of a Vector.
 */
Vector<float> Imag(const Vector<CS> &v) {
  return Imag((Matrix<CS>)v);
}

/**
 * Return the real part of a Vector.
 */
Vector<double> Real(const Vector<CD> &v) {
  return Real((Matrix<CD>)v);
}

/**
 * Return the imaginary part of a Vector.
 */
Vector<double> Imag(const Vector<CD> &v) {
  return Imag((Matrix<CD>)v);
}

/**
 * Join two vectors.
 */
template<typename T>
Vector<T> Join(const Vector<T> &vector1, const Vector<T> &vector2) {
  Vector<T> vector(vector1.size() + vector2.size());
  memcpy(vector.raw_ptr(), vector1.raw_ptr(), vector1.size() * sizeof(T));
  memcpy(vector.raw_ptr() + vector1.size(),
         vector2.raw_ptr(), vector2.size() * sizeof(T));
  vector.set_temporary(true);
  return vector;
}
}
#endif  // SRC_UTIL_VECTOR_BASIC_H_
