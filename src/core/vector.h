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

#ifndef SRC_CORE_VECTOR_H_
#define SRC_CORE_VECTOR_H_
#include <math.h>

#include <core/matrix.h>

#include <algorithm>
#include <vector>
#include <string>

namespace juzhen {

/////////////////////////////////////////////////////////////////////////////
/**
 *  Vector class
 */
template<typename DataType>
class Vector : public Matrix<DataType> {
 public:
  /**
   * Default constructor.
   */
  Vector();

  /**
   * Construct a Vector of size i.
   */
  Vector(size_t i);  // NOLINT

  /**
   * Construct a Vector from another vector.
   */
  template<typename T>
  Vector(const Vector<T> &m);  // NOLINT

  /**
   * Construct a Vector from another vector.
   */
  Vector(const Vector<DataType> &m);  // NOLINT

  /**
   * Construct a Vector from a Matrix's column-major raw array.
   */
  template<typename T>
  Vector(const Matrix<T> &m);  // NOLINT

  /**
   * Construct a Vector from a Matrix's column-major raw array.
   */
  Vector(const Matrix<DataType> &m);  // NOLINT

  /**
   * Construct a Vector from an array. The size of the Vector will be
   * the bigger of data's size and n.
   */
  template<typename T>
  Vector(const T *data, size_t n);

  /**
   * Construct a Vector from an existing Vector.
   */
  Vector(const std::vector<DataType> &v);  // NOLINT

  /**
   * Copy a Vector from another Vector.
   */
  Vector<DataType> &operator=(const Vector<DataType> &v);

  /**
   * Copy a Vector from a STL Vector.
   */
  Vector<DataType> &operator=(const std::vector<DataType> &v);

  /**
   * Return the vector itself.
   */
  Vector<DataType> &operator+();

  /**
   * Return the opposite vector.
   */
  Vector<DataType> operator-() const;

  /**
   * Assign all elements in the Vector to be rhs.
   */
  Vector<DataType> &Set(const DataType &rhs);

  /**
   * Resize a Vector.
   */
  void Resize(size_t n);

  /**
   * Resize a Vector to have nc*nr elements.
   */
  void Resize(size_t nc, size_t nr);

  /**
   * Return the size of the Vector.
   */
  inline size_t size() const;

  /**
   * Add a scalar and a vector.
   */
  Vector<DataType> operator+(const DataType &rhs) const;

  /**
   * Add two Vectors.
   */
  template<typename T>
  Vector<DataType> operator+(const Vector<T> &rhs) const;

  /**
   * Add a scalar.
   */
  Vector<DataType> &operator+=(const DataType &rhs);

  /**
   * Add two Vectors.
   */
  template<typename T>
  Vector<DataType> &operator+=(const Vector<T> &rhs);

  /**
   * Subtract a scalar.
   */
  Vector<DataType> operator-(const DataType &rhs) const;

  /**
   * Subtract two Vectors.
   */
  template<typename T>
  Vector<DataType> operator-(const Vector<T> &rhs) const;

  /**
   * Subtract a scalar.
   */
  Vector<DataType> &operator-=(const DataType &rhs);

   /**
   * Subtract two Vectors.
   */
  template<typename T>
  Vector<DataType> &operator-=(const Vector<T> &rhs);

  /**
   * Multiply a Vector with a real number.
   */
  Vector<DataType> operator*(const DataType &rhs) const;

  /**
   * Multiply a Vector with a real number.
   */
  Vector<DataType> &operator*=(const DataType &rhs);

  /**
   * Multiply a Vector and a Matrix. The Vector's size must be equal
   * to the Matrix's number of rows.
   */
  Vector<DataType> &operator*=(const Matrix<DataType> &rhs);

  /**
   * Divide a Vector by a constant.
   */
  Vector<DataType> operator/(const DataType &rhs) const;

  /**
   * Divide a Vector by a constant.
   */
  Vector<DataType> &operator/=(const DataType &rhs);

  /**
   * Get real part of a Vector.
   */
  template<typename T>
  Vector<T> Real() const;

  /**
   * Get imaginary part of a Vector.
   */
  template<typename T>
  Vector<T> Imag() const;

  /**
   * Get transpose of a Vector (just conjugate it).
   */
  Vector<DataType> Transpose() const;

  /**
   * Get hermitian of a Vector (just conjugate it).
   */
  Vector<DataType> Adjoint() const;

  /**
   * Get conjugate of a Vector.
   */
  Vector<DataType> Conjugate() const;

  /**
   * Swap two elements.
   */
  Vector<DataType> &Swap(size_t i1, size_t i2);

  /**
   * Return a sub-vector of [i1, i2)
   */
  Vector<DataType> Block(size_t i1, size_t i2);

  /**
   * Replace the sub-vector starting at i with another Vector v.
   */
  Vector<DataType> &Replace(size_t i, const Vector<DataType> &v);

  /**
   * Sort a Vector. This will give you error for a complex Vector.
   */
  Vector<DataType> &Sort();
};

typedef Vector<float> svector;
typedef Vector<double> dvector;
typedef Vector<CS> cvector;
typedef Vector<CD> zvector;

template<typename DataType>
Vector<DataType>::Vector() : Matrix<DataType>() { }

template<typename DataType>
Vector<DataType>::Vector(size_t i) : Matrix<DataType>(i, 1) { }

template<typename DataType>
template<typename T>
Vector<DataType>::Vector(const Vector<T> &v) {
  Matrix<DataType>::data_ptr_ = typename Matrix<DataType>::DataPtr(
      new DataArray<DataType>(v.raw_ptr(), v.size()));
  Matrix<DataType>::num_col_ = 1;
  Matrix<DataType>::num_row_ = v.size();
  Matrix<DataType>::raw_ptr_ = Matrix<DataType>::data_ptr_->data_ptr;
  Matrix<DataType>::temporary_ = v.temporary();
}

template<typename DataType>
Vector<DataType>::Vector(const Vector<DataType> &v) {
  if (v.temporary()) {
    Matrix<DataType>::data_ptr_ =
        (const_cast<Vector<DataType>&>(v)).data_ptr_;
    Matrix<DataType>::temporary_ = true;
  } else {
    Matrix<DataType>::data_ptr_ = typename Matrix<DataType>::DataPtr(
        new DataArray<DataType>(*(v.data_ptr_)));
    Matrix<DataType>::temporary_ = false;
  }
  Matrix<DataType>::num_col_ = 1;
  Matrix<DataType>::num_row_ = v.size();
  Matrix<DataType>::raw_ptr_ = Matrix<DataType>::data_ptr_->data_ptr;
}

template<typename DataType>
template<typename T>
Vector<DataType>::Vector(const Matrix<T> &m) {
  Matrix<DataType>::data_ptr_ = typename Matrix<DataType>::DataPtr(
      new DataArray<DataType>(m.raw_ptr(), m.size()));
  Matrix<DataType>::num_col_ = 1;
  Matrix<DataType>::num_row_ = m.size();
  Matrix<DataType>::raw_ptr_ = Matrix<DataType>::data_ptr_->data_ptr;
  Matrix<DataType>::temporary_ = m.temporary();
}

template<typename DataType>
Vector<DataType>::Vector(const Matrix<DataType> &m) {
  if (m.temporary()) {
    Matrix<DataType>::data_ptr_ =
        (const_cast<Matrix<DataType>&>(m)).data_ptr();
    Matrix<DataType>::temporary_ = true;
  } else {
    Matrix<DataType>::data_ptr_ = typename Matrix<DataType>::DataPtr(
      new DataArray<DataType>(m.raw_ptr(), m.size()));
    Matrix<DataType>::temporary_ = false;
  }
  Matrix<DataType>::num_col_ = 1;
  Matrix<DataType>::num_row_ = m.size();
  Matrix<DataType>::raw_ptr_ = Matrix<DataType>::data_ptr_->data_ptr;
}

template<typename DataType>
template<typename T>
Vector<DataType>::Vector(const T *data, size_t n)
  : Matrix<DataType>(data, n, 1) {}

template<typename DataType>
Vector<DataType>::Vector(const std::vector<DataType> &v)
  : Matrix<DataType>(&v[0], v.size(), 1) {
}

template<typename DataType>
Vector<DataType> &Vector<DataType>::operator=(const Vector<DataType> &v) {
  Matrix<DataType>::operator=(v);
  return *this;
}

template<typename DataType>
Vector<DataType> &Vector<DataType>::operator=(
    const std::vector<DataType> &v) {
  size_t s = v.size();
  Resize(s);
  DataType *p1;
  p1 = Matrix<DataType>::raw_ptr();
  for (size_t i = 0; i < s; i++) *(p1++) = v[i];
  return *this;
}

template<typename DataType>
Vector<DataType> &Vector<DataType>::operator+() {
  return *this;
}

template<typename DataType>
Vector<DataType> Vector<DataType>::operator-() const {
  return Matrix<DataType>::operator-();
}

template<typename DataType>
Vector<DataType> &Vector<DataType>::Set(const DataType &rhs) {
  Matrix<DataType>::Set(rhs);
  return *this;
}

template<typename DataType>
void Vector<DataType>::Resize(size_t n) {
  Matrix<DataType>::Resize(n, 1);
}

template<typename DataType>
void Vector<DataType>::Resize(size_t nc, size_t nr) {
  Matrix<DataType>::Resize(nc*nr, 1);
}

template<typename DataType>
inline size_t Vector<DataType>::size() const {
  return Matrix<DataType>::num_row_;
}

template<typename DataType>
Vector<DataType> Vector<DataType>::operator+(const DataType &rhs) const {
  return Matrix<DataType>::operator+(rhs);
}

template<typename DataType>
template<typename T>
Vector<DataType> Vector<DataType>::operator+(const Vector<T> &rhs) const {
  return Matrix<DataType>::operator+(rhs);
}

template<typename DataType>
Vector<DataType> &Vector<DataType>::operator+=(const DataType &rhs) {
  Matrix<DataType>::operator+=(rhs);
  return *this;
}

template<typename DataType>
template<typename T>
Vector<DataType> &Vector<DataType>::operator+=(const Vector<T> &rhs) {
  Matrix<DataType>::operator+=(rhs);
  return *this;
}

template<typename DataType>
Vector<DataType> Vector<DataType>::operator-(const DataType &rhs) const {
  return Matrix<DataType>::operator-(rhs);
}

template<typename DataType>
template<typename T>
Vector<DataType> Vector<DataType>::operator-(const Vector<T> &rhs) const {
  return Matrix<DataType>::operator-(rhs);
}

template<typename DataType>
Vector<DataType> &Vector<DataType>::operator-=(const DataType &rhs) {
  Matrix<DataType>::operator-=(rhs);
  return *this;
}

template<typename DataType>
template<typename T>
Vector<DataType> &Vector<DataType>::operator-=(const Vector<T> &rhs) {
  Matrix<DataType>::operator-=(rhs);
  return *this;
}

template<typename DataType>
Vector<DataType> Vector<DataType>::operator*(const DataType &rhs) const {
  return Matrix<DataType>::operator*(rhs);
}

template<typename DataType>
Vector<DataType> &Vector<DataType>::operator*=(const DataType &rhs) {
  Matrix<DataType>::operator*=(rhs);
  return *this;
}

template<typename DataType>
Vector<DataType> &Vector<DataType>::operator*=(const Matrix<DataType> &rhs) {
  Matrix<DataType> m = *this;
  (*this) = m.Transpose() * rhs;
  return *this;
}

template<typename DataType>
Vector<DataType> Vector<DataType>::operator/(const DataType &rhs) const {
  return Matrix<DataType>::operator/(rhs);
}

template<typename DataType>
Vector<DataType> &Vector<DataType>::operator/=(const DataType &rhs) {
  Matrix<DataType>::operator/=(rhs);
  return *this;
}

template<typename DataType>
Vector<DataType> Vector<DataType>::Transpose() const {
  return *this;
}

template<typename DataType>
Vector<DataType> Vector<DataType>::Adjoint() const {
  return Matrix<DataType>::Conjugate();
}

template<typename DataType>
Vector<DataType> Vector<DataType>::Conjugate() const {
  return Matrix<DataType>::Conjugate();
}

template<typename DataType>
Vector<DataType> &Vector<DataType>::Swap(size_t i1, size_t i2) {
  assert(i1 < size() && i2 < size());
  DataType *p = Matrix<DataType>::data_ptr_->data_ptr;
  DataType temp;
  temp = p[i1];
  p[i1] = p[i2];
  p[i2] = temp;
  return *this;
}

template<typename DataType>
Vector<DataType> Vector<DataType>::Block(size_t i1, size_t i2) {
  return Matrix<DataType>::Block(i1, 0, i2, 0);
}

template<typename DataType>
Vector<DataType> &Vector<DataType>::Replace(size_t i,
    const Vector<DataType> &v) {
  size_t s = size();
  for (size_t j = i, k = 0; j < s; j++, k++)
    (*this)[j] = v[k];
  return (*this);
}

template<typename DataType>
Vector<DataType> &Vector<DataType>::Sort() {
  std::vector<DataType> v(Matrix<DataType>::raw_ptr(),
                          Matrix<DataType>::raw_ptr()+size());
  std::sort(v.begin(), v.end());
  for (size_t i = 0; i < v.size(); i++) (*this)[i] = v[i];
  return (*this);
}
/////////////////////////////////////////////////////////////////////////////
}
#endif  // SRC_CORE_VECTOR_H_
