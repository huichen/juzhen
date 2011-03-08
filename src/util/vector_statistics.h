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

#ifndef SRC_UTIL_VECTOR_STATISTICS_H_
#define SRC_UTIL_VECTOR_STATISTICS_H_
#include <math.h>

#include <core/vector.h>

#include <algorithm>
#include <vector>
#include <string>

namespace mlcpp {
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
 * Return the average of a Vector's all elements.
 */
template<typename T>
T Average(const Vector<T> &vector) {
  T sum = 0;
  size_t endi = vector.size();
  if (endi == 0) return 0;
  for (size_t i = 0; i < endi; i++)
    sum += vector(i);
  return sum/endi;
}

/**
 * Find the Norm-2 of a Vector.
 */
double Norm(const Vector<CD> &v)  {
  double r = 0;
  CD *p = v.raw_ptr();
  size_t size = v.size();
  for (size_t i = 0; i < size; i++)
    r += abs2(*(p++));
  return sqrt(r);
}

float Norm(const Vector<CS> &v)  {
  float r = 0;
  CS *p = v.raw_ptr();
  size_t size = v.size();
  for (size_t i = 0; i < size; i++)
    r += abs2(*(p++));
  return sqrt(r);
}

double Norm(const Vector<double> &v)  {
  double r = 0;
  double *p = v.raw_ptr();
  size_t size = v.size();
  for (size_t i = 0; i < size; i++)
    r += abs2(*(p++));
  return sqrt(r);
}

float Norm(const Vector<float> &v)  {
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
  Vector<T> v2(v.size());
  std::vector<T> v1(v.raw_ptr(), v.raw_ptr()+v.size());
  std::sort(v1.begin(), v1.end());
  for (size_t i = 0; i < v1.size(); i++) v2[i] = v1[i];
  return v2;
}
}
#endif  // SRC_UTIL_VECTOR_STATISTICS_H_