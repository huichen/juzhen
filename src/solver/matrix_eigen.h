/*
+---------------------------------------------------------------------------+
|  Juzhen: C++ library for linear algebra                                           |
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

#ifndef SRC_SOLVER_MATRIX_EIGEN_H_
#define SRC_SOLVER_MATRIX_EIGEN_H_

#include <juzhen.h>

namespace juzhen {

/**
 * Solve eigen system and put eigen values into e, corresponding
 * left-eigen vectors into vl and corresponding right-eigen vectors
 * into vr. vl and vr will be square matrices.
 *
 * If you are linking to basic blas library only, this function is not
 * implemented yet.
 *
 * Set e to be empty matrix if it fails.
 */
template<typename DataType, typename T>
void EigenSolver(
    const Matrix<DataType> &mat,
    Matrix<Complex<T> > &e,
    Matrix<DataType> &vl,
    Matrix<DataType> &vr) {
  assert(mat.num_col() == mat.num_row());
  Matrix<DataType> m(mat);

  vl.Resize(mat.num_col(), mat.num_col());
  vr.Resize(mat.num_col(), mat.num_col());
  e.Resize(mat.num_col(), 1);

  char nl = &vl ? 'V' : 'N';
  char nr = &vr ? 'V' : 'N';

  if (geev<DataType>(nl, nr, mat.num_col(), m.raw_ptr(), mat.num_col(),
                 e.raw_ptr(), vl.raw_ptr(), mat.num_col(),
                 vr.raw_ptr(), mat.num_col()))
    e.Resize(0, 0);
}

/**
 * Solve eigen system and put eigen values into e, corresponding
 * right-eigen vectors into vr. vr will be square Matrix.
 *
 * If you are linking to basic blas library only, this function is not
 * implemented yet.
 *
 * Set e to be empty matrix if it fails.
 */
template<typename DataType, typename T>
void RightEigenSolver(
    const Matrix<DataType> &mat,
    Matrix<Complex<T> > &e,
    Matrix<DataType> &vr) {
  assert(mat.num_col() == mat.num_row());
  Matrix<DataType> m(mat);

  vr.Resize(mat.num_col(), mat.num_col());
  e.Resize(mat.num_col(), 1);

  char nl = 'N';
  char nr = 'V';

  if (geev<DataType>(nl, nr, mat.num_col(), m.raw_ptr(), mat.num_col(),
                 e.raw_ptr(), NULL, mat.num_col(),
                 vr.raw_ptr(), mat.num_col()))
    e.Resize(0, 0);
}

/**
 * Solve eigen system and put eigen values into e, corresponding
 * left-eigen vectors into vl. vl will be square Matrix.
 *
 * If you are linking to basic blas library only, this function is not
 * implemented yet.
 *
 * Set e to be empty matrix if it fails.
 */
template<typename DataType, typename T>
void LeftEigenSolver(
    const Matrix<DataType> &mat,
    Matrix<Complex<T> > &e,
    Matrix<DataType> &vl) {
  assert(mat.num_col() == mat.num_row());
  Matrix<DataType> m(mat);

  vl.Resize(mat.num_col(), mat.num_col());
  e.Resize(mat.num_col(), 1);

  char nl = 'V';
  char nr = 'N';

  if (geev<DataType>(nl, nr, mat.num_col(), m.raw_ptr(), mat.num_col(),
                 e.raw_ptr(), vl.raw_ptr(), mat.num_col(),
                 NULL, mat.num_col()))
    e.Resize(0, 0);
}
}
#endif  // SRC_SOLVER_MATRIX_EIGEN_H_
