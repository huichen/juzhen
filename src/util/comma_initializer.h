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

#ifndef SRC_UTIL_COMMA_INITIALIZER_H_
#define SRC_UTIL_COMMA_INITIALIZER_H_
#include <assert.h>

namespace juzhen {
/////////////////////////////////////////////////////////////////////////////
/**
 * Initialize a matrix with comma-separated values
 *
 * Example:
 *   dmatrix m(3,3);
 *   m << 1, 2, 3,
 *        4, 5, 6,
 *        7, 8, 9;
 * Matrix m is initialized row by row.
 */
template<typename T>
class MatrixCommaInitializer {
 public:
  MatrixCommaInitializer();
  MatrixCommaInitializer(const MatrixCommaInitializer<T> &init);  // NOLINT
  MatrixCommaInitializer(const Matrix<T> &matrix);  // NOLINT
  MatrixCommaInitializer<T> operator,(const T &value);  // NOLINT

 private:
  const Matrix<T> *super;
  size_t col;
  size_t row;
};

template<typename T>
MatrixCommaInitializer<T>::MatrixCommaInitializer()
    : super(NULL), col(0), row(0) {}

template<typename T>
MatrixCommaInitializer<T>::MatrixCommaInitializer(
    const MatrixCommaInitializer<T> &init)
    : super(init.super), col(init.col), row(init.row) {}

template<typename T>
MatrixCommaInitializer<T>::MatrixCommaInitializer(const Matrix<T> &matrix)
    : super(&matrix), col(0), row(0) {}

template<typename T>
MatrixCommaInitializer<T> MatrixCommaInitializer<T>::operator,(  // NOLINT
    const T &value) {
  const_cast<Matrix<T> *>(super)->operator()(row, col) = value;
  col++;
  if (col >= super->num_col()) {
    col = 0;
    row++;
  }

  if (row >= super->num_row()) {
    col = 0;
    row = 0;
  }
  return *this;
}
/////////////////////////////////////////////////////////////////////////////

template<typename T, typename T1>
MatrixCommaInitializer<T> operator<<(
    const Matrix<T> &matrix, const T1 &value) {
  return MatrixCommaInitializer<T>(matrix).operator,(value);  // NOLINT
}

/////////////////////////////////////////////////////////////////////////////
/**
 * Initialize a vector with comma-separated values
 *
 * Example:
 *   dmatrix m(6);
 *   m << 1, 2, 3, 4, 5, 6;
 */
template<typename T>
class VectorCommaInitializer {
 public:
  VectorCommaInitializer();
  VectorCommaInitializer(const VectorCommaInitializer<T> &init);  // NOLINT
  VectorCommaInitializer(const Vector<T> &matrix);  // NOLINT
  VectorCommaInitializer<T> operator,(const T &value);  // NOLINT

 private:
  const Vector<T> *super;
  size_t pos;
};

template<typename T>
VectorCommaInitializer<T>::VectorCommaInitializer()
    : super(NULL), pos(0) {}

template<typename T>
VectorCommaInitializer<T>::VectorCommaInitializer(
    const VectorCommaInitializer<T> &init)
    : super(init.super), pos(init.pos) {}

template<typename T>
VectorCommaInitializer<T>::VectorCommaInitializer(const Vector<T> &matrix)
    : super(&matrix), pos(0) {}

template<typename T>
VectorCommaInitializer<T> VectorCommaInitializer<T>::operator,(  // NOLINT
    const T &value) {
  const_cast<Vector<T> *>(super)->operator()(pos) = value;
  pos++;
  if (pos >= super->size()) {
    pos = 0;
  }
  return *this;
}
/////////////////////////////////////////////////////////////////////////////

template<typename T, typename T1>
VectorCommaInitializer<T> operator<<(
    const Vector<T> &matrix, const T1 &value) {
  return VectorCommaInitializer<T>(matrix).operator,(value);  // NOLINT
}
}
#endif  // SRC_UTIL_COMMA_INITIALIZER_H_
