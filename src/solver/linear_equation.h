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

#ifndef SRC_SOLVER_LINEAR_EQUATION_H_
#define SRC_SOLVER_LINEAR_EQUATION_H_

#include <mlcpp.h>

namespace mlcpp {

/**
 * Solve matrix equation A * X = B
 */
template<typename DataType>
Matrix<DataType> LinearSolver(
    const Matrix<DataType> &matrix_a,
    const Matrix<DataType> &matrix_b) {
  assert(matrix_a.num_col() == matrix_a.num_row());
  assert(matrix_b.num_col() == matrix_b.num_row());
  assert(matrix_a.num_col() == matrix_b.num_row());

  Matrix<DataType> mata;
  mata = matrix_a;

  Matrix<DataType> matrix_x;
  matrix_x = matrix_b;

  gesv<DataType>(mata.num_col(), mata.num_row(),
                 mata.raw_ptr(), mata.num_col(),
                 matrix_x.raw_ptr(), matrix_x.num_col());
  matrix_x.set_temporary(true);
  return matrix_x;
}
}
#endif  // SRC_SOLVER_LINEAR_EQUATION_H_
