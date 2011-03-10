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

#ifndef SRC_UTIL_VECTOR_STL_H_
#define SRC_UTIL_VECTOR_STL_H_
#include <math.h>

#include <core/vector.h>

#include <algorithm>
#include <vector>
#include <string>

namespace juzhen {

template<typename T>
std::ostream &operator<< (std::ostream &out, const Vector<T> &m) {
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
template<typename T>
std::string String(const Vector<T> &vector) {
  std::ostringstream out;
  out << vector;
  return out.str();
}

/**
 * Build a STL vector from a JUZHEN Vector.
 */
template<typename T>
std::vector<T> STLVector(const Vector<T> &vector) {
  std::vector<T> stl_vector(vector.raw_ptr(),
                          vector.raw_ptr()+vector.size());
  return stl_vector;
}
}
#endif  // SRC_UTIL_VECTOR_STL_H_
