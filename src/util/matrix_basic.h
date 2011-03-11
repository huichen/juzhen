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

#ifndef SRC_UTIL_MATRIX_BASIC_H_
#define SRC_UTIL_MATRIX_BASIC_H_
#include <core/matrix.h>

#include <sstream>
#include <string>

namespace juzhen {

/**
 * Return the Abs of a Matrix
 */
Matrix<float> Abs(const Matrix<float> &ma) {
  Matrix<float> m(ma);
  size_t endi = m.size();
  for (size_t i = 0; i < endi; i++)
    m(i) = Abs(m(i));
  m.set_temporary(true);
  return m;
}

Matrix<double> Abs(const Matrix<double> &ma) {
  Matrix<double> m(ma);
  size_t endi = m.size();
  for (size_t i = 0; i < endi; i++)
    m(i) = Abs(m(i));
  m.set_temporary(true);
  return m;
}

Matrix<float> Abs(const Matrix<CS> &ma) {
  Matrix<float> m(ma.num_row(), ma.num_col());
  size_t endi = ma.size();
  for (size_t i = 0; i < endi; i++)
    m(i) = Abs(ma(i));
  m.set_temporary(true);
  return m;
}

Matrix<double> Abs(const Matrix<CD> &ma) {
  Matrix<double> m(ma.num_row(), ma.num_col());
  size_t endi = ma.size();
  for (size_t i = 0; i < endi; i++)
    m(i) = Abs(ma(i));
  m.set_temporary(true);
  return m;
}

/**
 * Return the Abs-square of a Matrix
 */
Matrix<float> Abs2(const Matrix<float> &ma) {
  Matrix<float> m(ma);
  size_t endi = m.size();
  for (size_t i = 0; i < endi; i++)
    m(i) = Abs2(m(i));
  m.set_temporary(true);
  return m;
}

Matrix<double> Abs2(const Matrix<double> &ma) {
  Matrix<double> m(ma);
  size_t endi = m.size();
  for (size_t i = 0; i < endi; i++)
    m(i) = Abs2(m(i));
  m.set_temporary(true);
  return m;
}

Matrix<float> Abs2(const Matrix<CS> &ma) {
  Matrix<float> m(ma.num_row(), ma.num_col());
  size_t endi = ma.size();
  for (size_t i = 0; i < endi; i++)
    m(i) = Abs2(ma(i));
  m.set_temporary(true);
  return m;
}

Matrix<double> Abs2(const Matrix<CD> &ma) {
  Matrix<double> m(ma.num_row(), ma.num_col());
  size_t endi = ma.size();
  for (size_t i = 0; i < endi; i++)
    m(i) = Abs2(ma(i));
  m.set_temporary(true);
  return m;
}

/**
 * Return the real part of a Matrix
 */
Matrix<float> Real(const Matrix<float> &ma) {
  return ma;
}

/**
 * Return the imaginary part of a Matrix
 */
Matrix<float> Imag(const Matrix<float> &ma) {
  Matrix<float> m(ma);
  m.Clear();
  m.set_temporary(true);
  return m;
}

/**
 * Return the real part of a Matrix
 */
Matrix<double> Real(const Matrix<double> &ma) {
  return ma;
}

/**
 * Return the imaginary part of a Matrix
 */
Matrix<double> Imag(const Matrix<double> &ma) {
  Matrix<double> m(ma);
  m.Clear();
  m.set_temporary(true);
  return m;
}

/**
 * Return the real part of a Matrix
 */
Matrix<float> Real(const Matrix<CS> &ma) {
  Matrix<float> m(ma.num_row(), ma.num_col());
  size_t endi = ma.size();
  float *p1;
  CS *p2;
  p1 = m.raw_ptr();
  p2 = ma.raw_ptr();
  for (size_t i = 0; i < endi; i++)
    *(p1++) = (p2++)->real;
  m.set_temporary(true);
  return m;
}

/**
 * Return the imaginary part of a Matrix
 */
Matrix<float> Imag(const Matrix<CS> &ma) {
  Matrix<float> m(ma.num_row(), ma.num_col());
  size_t endi = ma.size();
  float *p1;
  CS *p2;
  p1 = m.raw_ptr();
  p2 = ma.raw_ptr();
  for (size_t i = 0; i < endi; i++)
    *(p1++) = (p2++)->imag;
  m.set_temporary(true);
  return m;
}

/**
 * Return the real part of a Matrix
 */
Matrix<double> Real(const Matrix<CD> &ma) {
  Matrix<double> m(ma.num_row(), ma.num_col());
  size_t endi = ma.size();
  double *p1;
  CD *p2;
  p1 = m.raw_ptr();
  p2 = ma.raw_ptr();
  for (size_t i = 0; i < endi; i++)
    *(p1++) = (p2++)->real;
  m.set_temporary(true);
  return m;
}

/**
 * Return the imaginary part of a Matrix
 */
Matrix<double> Imag(const Matrix<CD> &ma) {
  Matrix<double> m(ma.num_row(), ma.num_col());
  size_t endi = ma.size();
  double *p1;
  CD *p2;
  p1 = m.raw_ptr();
  p2 = ma.raw_ptr();
  for (size_t i = 0; i < endi; i++)
    *(p1++) = (p2++)->imag;
  m.set_temporary(true);
  return m;
}

/**
 * Return transpose of a Matrix
 */
template<typename T>
Matrix<T> Transpose(const Matrix<T> &m) {
  return m.Transpose();
}

/**
 * Return hermitian of a Matrix
 */
template<typename T>
Matrix<T> Adjoint(const Matrix<T> &m) {
  return m.Adjoint();
}

/**
 * Return conjugate of a Matrix
 */
template<typename T>
Matrix<T> Conjugate(const Matrix<T> &m) {
  return m.Conjugate();
}

template<typename T>
T Det(const Matrix<T> &m) {
  assert(m.num_col() == m.num_row());
  Matrix<T> matrix(m);
  return matrix_determinant(matrix.num_row(), matrix.raw_ptr());
}
}
#endif  // SRC_UTIL_MATRIX_BASIC_H_
