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

#ifndef SRC_UTIL_MATRIX_STL_H_
#define SRC_UTIL_MATRIX_STL_H_
#include <core/matrix.h>

#include <sstream>
#include <string>

namespace mlcpp {
/* stream operator overload */
template<typename DataType>
std::ostream &operator<< (std::ostream &out,
                          const Matrix<DataType> &m) {
  for (size_t i = 0; i < m.num_row(); i++) {
    for (size_t j = 0; j < m.num_col(); j++)
      out << m(i, j) << " ";
    if (i != m.num_row() - 1) out << std::endl;
  }
  return out;
}

/**
 * Get print form of a Matrix.
 */
template<typename DataType>
std::string String(const Matrix<DataType> &m) {
  std::ostringstream out;
  out << m;
  return out.str();
}
}
#endif  // SRC_UTIL_MATRIX_STL_H_
