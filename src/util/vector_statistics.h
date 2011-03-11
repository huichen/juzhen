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

#ifndef SRC_UTIL_VECTOR_STATISTICS_H_
#define SRC_UTIL_VECTOR_STATISTICS_H_
#include <math.h>

#include <core/vector.h>

#include <algorithm>
#include <vector>
#include <string>

namespace juzhen {
/**
 * Return the sum of a Vector's all elements.
 */
template<typename T>
T Sum(const Vector<T> &vector) {
  T sum = 0;
  size_t endi = vector.size();
  for (size_t i = 0; i < endi; i++)
    sum += vector(i);
  return sum;
}

/**
 * Return the product of a Vector's all elements.
 */
template<typename T>
T Prod(const Vector<T> &vector) {
  size_t endi = vector.size();
  if (endi == 0)
    return 0;
  T prod = 1;
  for (size_t i = 0; i < endi; i++)
    prod *= vector(i);
  return prod;
}

/**
 * Return the mean of a Vector's all elements.
 */
template<typename T>
T Mean(const Vector<T> &vector) {
  T sum = 0;
  size_t endi = vector.size();
  if (endi == 0) return 0;
  for (size_t i = 0; i < endi; i++)
    sum += vector(i);
  return sum/endi;
}

/**
 * Return the norm-square of a Vector.
 */
template<typename T>
double NormSquare(const Vector<T> &v) {
  double r = 0;
  size_t size = v.size();
  for (size_t i = 0; i < size; i++)
    r += Abs2(v(i));
  return r;
}

/**
 * Return the norm of a Vector.
 */
template<typename T>
double Norm(const Vector<T> &v) {
  return sqrt(NormSquare(v));
}

/**
 * Return the 1-norm of a Vector.
 */
template<typename T>
double NormOne(const Vector<T> &v) {
  double r = 0;
  size_t size = v.size();
  for (size_t i = 0; i < size; i++)
    r += Abs(v(i));
  return r;
}

/**
 * Return the infinity-norm of a Vector.
 */
template<typename T>
double NormInfinity(const Vector<T> &v) {
  size_t size = v.size();
  if (size == 0)
    return 0;
  double r = Abs(v(0));
  for (size_t i = 1; i < size; i++)
    if (r < Abs(v(i)))
      r = Abs(v(i));
  return r;
}

/**
 * Find the maximum element in a Vector.
 */
template<typename T>
T Max(const Vector<T> &vector) {
  size_t endi = vector.size();
  if (endi == 0)
    return 0;
  T max_value = vector(0);
  for (size_t i = 1; i < endi; i++)
    if (max_value < vector(i)) max_value = vector(i);
  return max_value;
}

/**
 * Find the minimum element in a Vector.
 */
template<typename T>
T Min(const Vector<T> &vector) {
  size_t endi = vector.size();
  if (endi == 0)
    return 0;
  T min_value = vector(0);
  for (size_t i = 1; i < endi; i++)
    if (min_value > vector(i)) min_value = vector(i);
  return min_value;
}

/**
 * Sort a Vector.
 */
template<typename T>
Vector<T> Sort(const Vector<T> &v) {
  Vector<T> v1(v);
  std::vector<T> v2(v1.raw_ptr(), v1.raw_ptr()+v1.size());
  std::sort(v2.begin(), v2.end());
  for (size_t i = 0; i < v1.size(); i++) v1[i] = v2[i];
  v1.set_temporary(true);
  return v1;
}

/**
 * Return covariance of two vectors.
 */
template<typename T>
T Cov(const Vector<T> &vector1, const Vector<T> &vector2) {
  assert(vector1.size() == vector2.size());
  size_t endi = vector1.size();
  if (endi == 0)
    return 0;
  T mean1 = Mean(vector1);
  T mean2 = Mean(vector2);
  T covar = 0;
  for (size_t i = 0; i < endi; i++)
    covar += (vector1(i) - mean2) * (vector2(i) - mean2);
  return covar/endi;
}

/**
 * Return variance of a vector.
 */
template<typename T>
T Var(const Vector<T> &vector) {
  return Cov(vector, vector);
}

/**
 * Return standard deviation of a vector.
 */
template<typename T>
T StdDev(const Vector<T> &vector) {
  assert(0);
}

template<>
float StdDev(const Vector<float> &vector) {
  return sqrt(Var(vector));
}

template<>
double StdDev(const Vector<double> &vector) {
  return sqrt(Var(vector));
}

/**
 * Return correlation coefficient of two vectors.
 */
template<typename T>
T CorrCoeff(const Vector<T> &vector1, const Vector<T> &vector2) {
  return Cov(vector1, vector2) / (StdDev(vector1) * StdDev(vector2));
}
}
#endif  // SRC_UTIL_VECTOR_STATISTICS_H_
