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

#ifndef SRC_SOLVER_MATRIX_INVERSE_H_
#define SRC_SOLVER_MATRIX_INVERSE_H_

#include <juzhen.h>

namespace juzhen {

/**
 * Inverse a matrix
 * X = A^(-1)
 *
 * Return empty matrix if it fails.
 */
template<typename DataType>
Matrix<DataType> Inverse(
    const Matrix<DataType> &matrix_a) {
  assert(matrix_a.num_col() == matrix_a.num_row());

  Matrix<DataType> matrix_x;
  matrix_x = matrix_a;

  if (matrix_inverse<DataType>(matrix_x.num_col(), matrix_x.num_row(),
                 matrix_x.raw_ptr(), matrix_x.num_col()))
    return Matrix<DataType>();
  matrix_x.set_temporary(true);
  return matrix_x;
}
}
#endif  // SRC_SOLVER_MATRIX_INVERSE_H_
