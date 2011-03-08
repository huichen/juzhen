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

#ifndef SRC_UTIL_VECTOR_MATH_H_
#define SRC_UTIL_VECTOR_MATH_H_
#include <math.h>

#include <core/vector.h>

#include <algorithm>
#include <vector>
#include <string>

namespace mlcpp {
/**
 * Multiply a real number and a Vector.
 */
template<typename DataType>
Vector<DataType> operator*(double lhs, const Vector<DataType> &ma) {
  if (ma.temporary_) {
    DataType *p = ma.raw_ptr();
    size_t maxi = ma.size();
    for (size_t i = 0; i < maxi; i++)
      *(p++) *= lhs;
    return ma;
  } else {
    Vector<DataType> m(ma.size());
    DataType *p1 = ma.raw_ptr();
    DataType *p2 = m.raw_ptr();
    size_t maxi = ma.size();
    for (size_t i = 0; i < maxi; i++)
      *(p2++) = *(p1++)*lhs;
    return m;
  }
}

/**
 * Multiply a complex number and a Vector.
 */
template<typename DataType>
Vector<DataType> operator*(const CD &lhs, const Vector<DataType> &ma) {
  if (ma.temporary_) {
    DataType *p = ma.raw_ptr();
    size_t maxi = ma.size();
    for (size_t i = 0; i < maxi; i++)
      *(p++) *= lhs;
    return ma;
  } else {
    Vector<DataType> m(ma.size());
    DataType *p1 = ma.raw_ptr();
    DataType *p2 = m.raw_ptr();
    size_t maxi = ma.size();
    for (size_t i = 0; i < maxi; i++)
      *(p2++) = *(p1++)*lhs;
    return m;
  }
}

/**
 * Multiply a Matrix and a Vector.
 */
template<typename DataType>
Vector<DataType> operator*(
    const Matrix<DataType> &ma,
    const Vector<DataType> &v) {
  return ma * ((Matrix<DataType>&)v);
}

/**
 * Multiply a Matrix and a Vector.
 */
template<typename DataType>
Vector<DataType> operator*(
    const Vector<DataType> &v,
    const Matrix<DataType> &ma) {
  return ma * ((Matrix<DataType>&)v);
}

/**
 * Multiply a Matrix and a Vector.
 */
template<typename DataType>
Vector<DataType> &operator*=(
    Matrix<DataType> &ma,
    const Vector<DataType> &v) {
  ma *= (Matrix<DataType>&)v;
  return ma;
}

/**
 * Vector dot product.
 */
template<typename DataType>
DataType operator*(const Vector<DataType> &v1, const Vector<DataType> &v2) {
  assert(v1.size() == v2.size());
  DataType res = 0;
  size_t s = v1.size();
  DataType *p1 = v1.raw_ptr();
  DataType *p2 = v2.raw_ptr();
  for (size_t i = 0; i < s; i++) res+=*(p1++)*(*(p2++));
  return res;
}
}
#endif  // SRC_UTIL_VECTOR_MATH_H_
