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

#ifndef SRC_UTIL_VECTOR_MATH_H_
#define SRC_UTIL_VECTOR_MATH_H_
#include <math.h>

#include <core/vector.h>

#include <algorithm>
#include <vector>
#include <string>

namespace juzhen {

template<typename T>
Vector<T> operator+(double lhs, const Vector<T> &v) {
  return v + lhs;
}

template<typename T>
Vector<T> operator+(const CD &lhs, const Vector<T> &v) {
  return v + lhs;
}

template<typename T>
Vector<T> operator-(double lhs, const Vector<T> &v) {
  return (-v) + lhs;
}

template<typename T>
Vector<T> operator-(const CD &lhs, const Vector<T> &v) {
  return (-v) + lhs;
}

/**
 * Multiply a real number and a Vector.
 */
template<typename T>
Vector<T> operator*(double lhs, const Vector<T> &ma) {
  if (ma.temporary()) {
    T *p = ma.raw_ptr();
    size_t maxi = ma.size();
    for (size_t i = 0; i < maxi; i++)
      *(p++) *= lhs;
    return ma;
  } else {
    Vector<T> m(ma.size());
    T *p1 = ma.raw_ptr();
    T *p2 = m.raw_ptr();
    size_t maxi = ma.size();
    for (size_t i = 0; i < maxi; i++)
      *(p2++) = *(p1++)*lhs;
    return m;
  }
}

/**
 * Multiply a complex number and a Vector.
 */
template<typename T>
Vector<T> operator*(const CD &lhs, const Vector<T> &ma) {
  if (ma.temporary()) {
    T *p = ma.raw_ptr();
    size_t maxi = ma.size();
    for (size_t i = 0; i < maxi; i++)
      *(p++) *= lhs;
    return ma;
  } else {
    Vector<T> m(ma.size());
    T *p1 = ma.raw_ptr();
    T *p2 = m.raw_ptr();
    size_t maxi = ma.size();
    for (size_t i = 0; i < maxi; i++)
      *(p2++) = *(p1++)*lhs;
    return m;
  }
}

/**
 * Multiply a Matrix and a Vector.
 */
template<typename T>
Vector<T> operator*(
    const Matrix<T> &ma,
    const Vector<T> &v) {
  return ma * ((Matrix<T>&)v);
}

/**
 * Multiply a Matrix and a Vector.
 */
template<typename T>
Vector<T> operator*(
    const Vector<T> &v,
    const Matrix<T> &ma) {
  return ma * ((Matrix<T>&)v);
}

/**
 * Multiply a Matrix and a Vector.
 */
template<typename T>
Vector<T> &operator*=(
    Matrix<T> &ma,
    const Vector<T> &v) {
  ma *= (Matrix<T>&)v;
  return ma;
}

/**
 * Vector dot (inner) product.
 */
template<typename T>
T operator*(const Vector<T> &v1, const Vector<T> &v2) {
  assert(v1.size() == v2.size());
  T res = 0;
  size_t s = v1.size();
  T *p1 = v1.raw_ptr();
  T *p2 = v2.raw_ptr();
  for (size_t i = 0; i < s; i++) res+=*(p1++)*Conjugate(*(p2++));
  return res;
}

/**
 * Vector cross product.
 */
template<typename T>
Matrix<T> OuterProduct(const Vector<T> &vector1, const Vector<T> &vector2) {
  Matrix<T> matrix(vector1.size(), vector2.size());
  size_t endi = vector1.size();
  size_t endj = vector2.size();
  for (size_t i = 0; i < endi; i++)
    for (size_t j = 0; j < endj; j++)
      matrix(i, j) = vector1(i) * vector2(j);
  matrix.set_temporary(true);
  return matrix;
}

/**
 * Vector cross product.
 */
template<typename T>
Vector<T> CrossProduct(const Vector<T> &vector1, const Vector<T> &vector2) {
  assert(vector1.size() == 3);
  assert(vector2.size() == 3);
  Vector<T> vector(3);
  vector(0) = vector1(1) * vector2(2) - vector1(2) * vector2(1);
  vector(1) = vector1(2) * vector2(0) - vector1(0) * vector2(2);
  vector(2) = vector1(0) * vector2(1) - vector1(1) * vector2(0);
  vector.set_temporary(true);
  return vector;
}
}
#endif  // SRC_UTIL_VECTOR_MATH_H_
