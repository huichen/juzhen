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

#ifndef SRC_UTIL_VECTOR_UTIL_H_
#define SRC_UTIL_VECTOR_UTIL_H_
#include <math.h>

#include <core/vector.h>

#include <algorithm>
#include <vector>
#include <string>

namespace mlcpp {

/**
 * Return the transpose of a Vector.
 */
template<typename DataType>
Vector<DataType> Transpose(const Vector<DataType> &v) {
  return v;
}

/**
 * Return the hermitian of a Vector.
 */
template<typename DataType>
Vector<DataType> Adjoint(const Vector<DataType> &v) {
  return v.Adjoint();
}

/**
 * Return the conjugate of a Vector.
 */
template<typename DataType>
Vector<DataType> Conjugate(const Vector<DataType> &v) {
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
 * Find the Sum of a Vector's all elements.
 */
template<typename DataType>
DataType Sum(const Vector<DataType> &v) {
  return v.Sum();
}

/**
 * Find the norm-2 of a Vector.
 */
double norm(const Vector<CD> &v)  {
  double r = 0;
  CD *p = v.raw_ptr();
  size_t size = v.size();
  for (size_t i = 0; i < size; i++)
    r += abs2(*(p++));
  return sqrt(r);
}

float norm(const Vector<CS> &v)  {
  float r = 0;
  CS *p = v.raw_ptr();
  size_t size = v.size();
  for (size_t i = 0; i < size; i++)
    r += abs2(*(p++));
  return sqrt(r);
}

double norm(const Vector<double> &v)  {
  double r = 0;
  double *p = v.raw_ptr();
  size_t size = v.size();
  for (size_t i = 0; i < size; i++)
    r += abs2(*(p++));
  return sqrt(r);
}

float norm(const Vector<float> &v)  {
  double r = 0;
  float *p = v.raw_ptr();
  size_t size = v.size();
  for (size_t i = 0; i < size; i++)
    r += abs2(*(p++));
  return sqrt(r);
}

/**
 * Find the maximum element in a Vector.
 */
template<typename DataType>
DataType Max(const Vector<DataType> &v) {
  return v.Max();
}

/**
 * Sort a Vector.
 */
template<typename DataType>
Vector<DataType> Sort(const Vector<DataType> &v) {
  Vector<DataType> v2(v.size());
  std::vector<DataType> v1(v.raw_ptr(), v.raw_ptr()+v.size());
  std::sort(v1.begin(), v1.end());
  for (size_t i = 0; i < v1.size(); i++) v2[i] = v1[i];
  return v2;
}

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

/**
 * Multiply a Matrix and a Vector.
 */
template<typename DataType>
Vector<DataType> operator*(
    const Matrix<DataType> &ma,
    const Vector<DataType> &v) {
  return ma*((Matrix<DataType>&)v);
}

/**
 * Multiply a Matrix and a Vector.
 */
template<typename DataType>
Vector<DataType> operator*(
    const Vector<DataType> &v,
    const Matrix<DataType> &ma) {
  return ma*((Matrix<DataType>&)v);
}

/**
 * Multiply a Matrix and a Vector.
 */
template<typename DataType>
Vector<DataType> &operator*=(
    Matrix<DataType> &ma,
    const Vector<DataType> &v) {
  ma *= (Matrix<DataType>&) v;
  return ma;
}

template<typename DataType>
std::ostream &operator<< (std::ostream &out, const Vector<DataType> &m) {
  out << "{";
  for (size_t i = 0; i < m.size(); i++) {
    if (i != m.size()-1)
      out << m(i) << ", ";
    else
      out << m(i);
  }
  out << "}";
  return out;
}

/**
 * Return string form of a Matrix.
 */

template<typename DataType>
std::string OutputToString(const Vector<DataType> &m) {
  std::ostringstream out;
  out << m;
  return out.str();
}
}
/////////////////////////////////////////////////////////////////////////////
#endif  // SRC_UTIL_VECTOR_UTIL_H_
