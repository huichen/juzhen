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

#ifndef SRC_CORE_IDENTITY_MATRIX_H_
#define SRC_CORE_IDENTITY_MATRIX_H_
#include <core/matrix.h>

namespace juzhen {

/////////////////////////////////////////////////////////////////////////////
/**
 * Identity Matrix
 */
template<typename DataType>
class IdentityMatrix : public Matrix<DataType> {
 public:
  /**
   * Construct an identity Matrix of n x n.
   */
  IdentityMatrix(size_t n);  // NOLINT
};

typedef IdentityMatrix<float> identity_smatrix;
typedef IdentityMatrix<double> identity_dmatrix;
typedef IdentityMatrix<CS> identity_cmatrix;
typedef IdentityMatrix<CD> identity_zmatrix;

template<typename DataType>
IdentityMatrix<DataType>::IdentityMatrix(size_t n) {
  Matrix<DataType>::data_ptr_ = typename Matrix<DataType>::DataPtr(
                                           new DataArray<DataType>(n*n));
  Matrix<DataType>::num_col_ = n;
  Matrix<DataType>::num_row_ = n;
  Matrix<DataType>::raw_ptr_ = Matrix<DataType>::data_ptr_->data_ptr;

  size_t endi = n*n;
  DataType *p = Matrix<DataType>::raw_ptr();
  for (size_t i = 0; i < endi; i++)
    *(p++) = 0;

  p = Matrix<DataType>::raw_ptr();
  for (size_t i = 0; i < n; i++)
    *(p+i*(n+1)) = 1.;
}
/////////////////////////////////////////////////////////////////////////////
}
#endif  // SRC_CORE_IDENTITY_MATRIX_H_
