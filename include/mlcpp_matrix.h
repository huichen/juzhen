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

#ifndef MLCPP_MATRIX_H_  // NOLINT
#define MLCPP_MATRIX_H_
#include <sstream>
#include <string>

#include "mlcpp_complex.h"
#include "mlcpp_dataarray.h"
#include "mlcpp_adaptor.h"

namespace mlcpp {

/////////////////////////////////////////////////////////////////////////////
/**
 * The basic Matrix class that all variant Matrix classes are derived from
 */
template<typename DataType>
class Matrix {
 public:
  /**
   * Deconstructor.
   */
  ~Matrix() {}

  /**
   * Default constructor.
   */
  Matrix();

  /**
   * Construct an Matrix of nr rows and nc columns. Memory is allocated but
   * Matrix's items have undetermined values.
   */
  Matrix(size_t nr, size_t nc);

  /**
   * Construct an Matrix of nr rows and nc columns from a raw array. nr*nc
   * numbers will be allocated in memory. If data's size is larger than nr*nc
   * only first nr*nc numbers are used in column-major order. If data's
   * size is smaller than nr*nc, all numbers from data will be used and the
   * Matrix's rests number have undetermined values.
   */
  template<typename T>
  Matrix(const T *data, size_t nr, size_t nc);

  /**
   * Construct from another Matrix. All numbers are copied.
   */
  template<typename T>
  Matrix(const Matrix<T> &m);  // NOLINT

  /**
   * Construct from another Matrix. All numbers are copied.
   */
  explicit Matrix(const Matrix<DataType> &m);

  /**
   * Copy from another Matrix. New memory is allocated if it's necessary.
   */
  Matrix<DataType>& operator= (const Matrix<DataType> &rhs);

  /**
   * Assign all numbers in the Matrix to be rhs.
   */
  Matrix<DataType>& operator= (const DataType &rhs);

  /**
   * Check if two matrices are equal. The two matrices must have the same
   * numbers of rows and columns if it returns true.
   */
  bool operator== (const Matrix<DataType> &rhs);

  /**
   * Check if two matrices are not equal. Returns false if the two matrices
   * have different numbers of rows or columns.
   */
  bool operator!= (const Matrix<DataType> &rhs);

  /***********************************************
   *  interface to private data
   **********************************************/

  /**
   * Return number of rows that the Matrix has.
   */
  size_t num_row() const;

  /**
   * Return number of columns that the Matrix has.
   */
  size_t num_col() const;

  /**
   * Returns the pointer of raw array. This is pretty useful when calling
   * low level blas/lapack functions.
   */
  inline DataType * raw_ptr() const;

  /**
   * Get a Matrix's temporary flag.
   */
  inline bool temporary() const;

  /**
   * Set a Matrix's temporary flag.
   */
  inline void set_temporary(bool t);

  /**
   * Resizing or reshaping the Matrix. If the resized size is larger than
   * the original, new memory will be allocated and data will be copied to
   * the new position.
   */
  void Resize(size_t nr, size_t nc);

  /**
   * Set all elements to be zero.
   */
  inline void Clear();

  /**
   * Find the reference of item located at ith row and jth column.
   */
  inline DataType& operator()(size_t i, size_t j);

  /**
   * Find the reference of item located at ith row and jth column.
   */
  inline DataType operator()(size_t i, size_t j) const;

  /**
   * Find the reference of item located at ith position in the raw data array.
   * The array is saved in a column-major order in memory.
   */
  inline DataType& operator()(size_t i);

  /**
   * Find the reference of item located at ith position in the raw data array.
   * The array is saved in a column-major order in memory.
   */
  inline DataType& operator()(size_t i) const;

  /**
   * Find the reference of item located at ith position in the raw data array.
   * The array is saved in a column-major order in memory.
   */
  inline DataType& operator[](size_t i);

  /**
   * Find the reference of item located at ith position in the raw data array.
   * The array is saved in a column-major order in memory.
   */
  inline DataType& operator[](size_t i) const;

  /****************************
   * Arithmetic operations
   ***************************/

  /**
   * Add two matrices
   */
  template<typename T>
  Matrix<DataType> operator+(const Matrix<T>& rhs) const;

  /**
   * Add two matrices
   */
  Matrix<DataType> operator+(const Matrix<DataType>& rhs) const;

  /**
   * Add one matrices upto another Matrix.
   */
  template<typename T>
  Matrix<DataType>& operator+=(const Matrix<T>& rhs);

  /**
   * Subtract one Matrix and another Matrix.
   */
  template<typename T>
  Matrix<DataType> operator-(const Matrix<T>& rhs) const;

  /**
   * Subtract one Matrix and another Matrix.
   */
  Matrix<DataType> operator-(const Matrix<DataType>& rhs) const;

  /**
   * Subtract one Matrix from another Matrix.
   */
  template<typename T>
  Matrix<DataType>& operator-=(const Matrix<T>& rhs);

  /**
   * Multiply each element in the Matrix by a real number and
   * get a new Matrix.
   */
  Matrix<DataType> operator*(double rhs) const;

  /**
   * Multiply each element in the Matrix by a complex number and
   * get a new Matrix.
   */
  Matrix<DataType> operator*(const CD &rhs) const;

  /**
   * Multiply each element in the Matrix by a real number.
   */
  Matrix<DataType> & operator*=(double rhs);

  /**
   * Multiply each element in the Matrix by a complex number.
   */
  Matrix<DataType> & operator*=(const CD &rhs);

  /**
   * Multiply two matrices and get a new Matrix. Number of columns of the first
   * Matrix must equals to the number of columns of the second Matrix.
   */
  Matrix<DataType> operator*(const Matrix<DataType>& rhs) const;

  /**
   * Multiply the Matrix with a new Matrix. Number of columns of the first
   * Matrix must equals to the number of columns of the second Matrix.
   */
  Matrix<DataType>& operator*=(const Matrix<DataType>& rhs);

  /**
   * Divide each element in the Matrix by a constant to get a new Matrix.
   */
  template<typename T>
  Matrix<DataType> operator/(const T &rhs) const;

  /**
   * Divide each element in the Matrix by a constant.
   */
  template<typename T>
  Matrix<DataType>& operator/=(const T &rhs);

  /*************************************
   *  Matrix specific operations
   ************************************/
  /**
   * Get the Matrix's transposition.
   */
  Matrix<DataType> Transpose()const;

  /**
   * Get the Matrix's hermitian.
   */
  Matrix<DataType> Adjoint() const;

  /**
   * Get the Matrix's conjugate.
   */
  Matrix<DataType> Conjugate()const;

  /**
   * Get a sub-block of the Matrix between r1th (containing) row and r2th
   * (not containing) row, and between c1th (containing) column and c2th
   * (not containing) column.
   */
  Matrix<DataType> Block(size_t r1, size_t r2, size_t c1, size_t c2) const;

  /**
   * Replace a sub-block of the Matrix with another Matrix value. The
   * sub-block starts at (r-th row, c-th column) at its upper-left corner,
   * and is replace by values from mt's (0,0) position.
   */
  Matrix<DataType> &Replace(size_t r, size_t c, const Matrix<DataType> &mt);

  /**
   * Return the Matrix's c-th column.
   */
  inline Matrix<DataType> GetCol(size_t c) const;

  /**
   * Return the Matrix's r-th row.
   */
  Matrix<DataType> GetRow(size_t r) const;

  /**
   * Swap two cols.
   */
  Matrix<DataType> &SwapCol(size_t c1, size_t c2);

  /**
   * Swap two rows.
   */
  Matrix<DataType> &SwapRow(size_t r1, size_t r2);

  /**
   * Solve eigen system and put eigen values into e, corresponding
   * left-eigen vectors into vl and corresponding right-eigen vectors
   * into vr. vl and vr will be square matrices.
   *
   * If you are linking to basic blas library only, this function is not
   * implemented yet.
   */
  template<typename T>
  void EigenSolver(
      Matrix<Complex<T> > &e,
      Matrix<DataType> &vl,
      Matrix<DataType> &vr) const;

  /**
   * Solve eigen system and put eigen values into e, corresponding
   * left-eigen vectors into vl. vl will be square Matrix.
   *
   * If you are linking to basic blas library only, this function is not
   * implemented yet.
   */
  template<typename T>
  void RightEigenSolver(
      Matrix<Complex<T> > &e,
      Matrix<DataType> &vr) const;

  /**
   * Solve eigen system and put eigen values into e, corresponding
   * right-eigen vectors into vr. vr will be square Matrix.
   *
   * If you are linking to basic blas library only, this function is not
   * implemented yet.
   */
  template<typename T>
  void LeftEigenSolver(
      Matrix<Complex<T> > &e,
      Matrix<DataType> &vl) const;

  typedef std::auto_ptr<DataArray<DataType> > DataPtr;

 protected:
  /**
   * Pointer to raw array.
   */
  DataType * raw_ptr_;

  /**
   * Auto pointer to raw_ptr_.
   */
  DataPtr data_ptr_;

  /**
   * Number of columns.
   */
  size_t num_col_;

  /**
   * Number of rows.
   */
  size_t num_row_;

  /**
   * If the Matrix is temporary.
   */
  bool temporary_;
};

typedef Matrix<float> smatrix;
typedef Matrix<double> dmatrix;
typedef Matrix<CS> cmatrix;
typedef Matrix<CD> zmatrix;

template<typename DataType>
Matrix<DataType>::Matrix() : num_col_(0), num_row_(0) {
  temporary_ = false;
  raw_ptr_ = NULL;
}

template<typename DataType>
Matrix<DataType>::Matrix(size_t nr, size_t nc) : num_col_(nc), num_row_(nr) {
  data_ptr_ = DataPtr(new DataArray<DataType>(num_col_*num_row_));
  temporary_ = false;
  raw_ptr_ = data_ptr_->data_ptr;
}

template<typename DataType>
template<typename T>
Matrix<DataType>::Matrix(const T *data, size_t nr, size_t nc)
  : num_col_(nc), num_row_(nr) {
  data_ptr_ =DataPtr(new DataArray<DataType>(data, nr*nc));
  temporary_ = false;
  raw_ptr_ = data_ptr_->data_ptr;
}

template<typename DataType>
template<typename T>
Matrix<DataType>::Matrix(const Matrix<T> &m)
  : num_col_(m.num_col()), num_row_(m.num_row()) {
  data_ptr_ =DataPtr(
      new DataArray<DataType>(m.raw_ptr(), m.num_row()*m.num_col()));
  temporary_ = false;
  raw_ptr_ = data_ptr_->data_ptr;
}

template<typename DataType>
Matrix<DataType>::Matrix(const Matrix<DataType> &m)
  : num_col_(m.num_col()), num_row_(m.num_row()) {
  if (m.temporary_) {
    data_ptr_ = (const_cast<Matrix<DataType>& >(m)).data_ptr_;
    temporary_ = true;
  } else {
    data_ptr_ =DataPtr(new DataArray<DataType>(*(m.data_ptr_)));
    temporary_ = false;
  }
  raw_ptr_ = data_ptr_->data_ptr;
}

template<typename DataType>
Matrix<DataType>&
Matrix<DataType>::operator=(const Matrix<DataType> &rhs) {
  if (&rhs == this) return *this;
  if (rhs.temporary_) {
    data_ptr_ = (const_cast<Matrix<DataType>& >(rhs)).data_ptr_;
    num_col_ = rhs.num_col();
    num_row_ = rhs.num_row();
    raw_ptr_ = data_ptr_->data_ptr;
    return *this;
  }
  if ((raw_ptr_ == NULL) || data_ptr_->size < rhs.num_col()*rhs.num_row())
    data_ptr_ = DataPtr(new DataArray<DataType>(*(rhs.data_ptr_)));
  else
    memcpy(data_ptr_->data_ptr,
           rhs.raw_ptr(),
           rhs.num_col()*rhs.num_row()*sizeof(DataType));
  num_col_ = rhs.num_col();
  num_row_ = rhs.num_row();
  raw_ptr_ = data_ptr_->data_ptr;
  return *this;
}

template<typename DataType>
Matrix<DataType>&
Matrix<DataType>::operator=(const DataType &rhs) {
  size_t endi = num_col_*num_row_;
  DataType *p = data_ptr_->data_ptr;
  for (size_t i = 0; i < endi; i++) *(p++)=rhs;
  return *this;
}

template<typename DataType>
bool Matrix<DataType>::operator==(const Matrix<DataType> &rhs) {
  if (&rhs == this) return true;
  if (num_row_ != rhs.num_row() || num_col_ != rhs.num_col()) return false;
  size_t endi = num_col_*num_row_;
  DataType *p1 = data_ptr_->data_ptr;
  DataType *p2 = rhs.raw_ptr();
  for (size_t i = 0; i < endi; i++)
    if (*(p1++) != *(p2++)) return false;
  return true;
}

template<typename DataType>
bool Matrix<DataType>::operator!=(const Matrix<DataType> &rhs) {
  return !(operator == (rhs));
}

/* interface to private data */
template<typename DataType>
size_t Matrix<DataType>::num_row() const {
  return num_row_;
}

template<typename DataType>
size_t Matrix<DataType>::num_col() const {
  return num_col_;
}

template<typename DataType>
inline DataType * Matrix<DataType>::raw_ptr() const {
//  return data_ptr_->data_ptr;
  return raw_ptr_;
}

/* resizing/reshaping Matrix */
template<typename DataType>
void Matrix<DataType>::Resize(size_t nr, size_t nc) {
  if (!raw_ptr_ || data_ptr_->size < nc * nr) {
    DataPtr
      newData(new DataArray<DataType>(*data_ptr_, nr*nc));
    data_ptr_ = newData;
    raw_ptr_ = data_ptr_->data_ptr;
  }
  num_col_ = nc;
  num_row_ = nr;
}

template<typename DataType>
inline bool Matrix<DataType>::temporary() const {
  return temporary_;
}

template<typename DataType>
inline void Matrix<DataType>::set_temporary(bool t) {
  temporary_ = t;
};

template<typename DataType>
inline void Matrix<DataType>::Clear() {
  if (data_ptr_->data_ptr)
    memset(data_ptr_->data_ptr, 0, data_ptr_->size*sizeof(DataType));
}

/* operator overloading */
template<typename DataType>
inline DataType& Matrix<DataType>::operator()(size_t i, size_t j) {
  assert(i < num_row_ && j < num_col_);
  return raw_ptr_[j*num_row_ + i];
}

template<typename DataType>
inline DataType Matrix<DataType>::operator()(size_t i, size_t j) const {
  assert(i < num_row_ && j < num_col_);
  return raw_ptr_[j*num_row_ + i];
}

template<typename DataType>
inline DataType& Matrix<DataType>::operator()(size_t i) {
  return raw_ptr_[i];
}

template<typename DataType>
inline DataType& Matrix<DataType>::operator()(size_t i) const {
  return raw_ptr_[i];
}

template<typename DataType>
inline DataType& Matrix<DataType>::operator[](size_t i) {
  return raw_ptr_[i];
}

template<typename DataType>
inline DataType& Matrix<DataType>::operator[](size_t i) const {
  return raw_ptr_[i];
}

// arithmetic

template<typename DataType>
template<typename T>
Matrix<DataType> Matrix<DataType>::operator+(const Matrix<T>& rhs) const {
  assert(num_col_ == rhs.num_col() && num_row_ == rhs.num_row());
  if (temporary_) {
    (const_cast<Matrix<DataType>* >(this))->operator+=(rhs);
    return *this;
  } else {
    Matrix<DataType> m(num_row_, num_col_);
    size_t endi = num_row_*num_col_;
    DataType *p1, *p2, *p3;
    p1 = m.raw_ptr();
    p2 = raw_ptr();
    p3 = rhs.raw_ptr();
    for (size_t i = 0; i < endi; i++)
      *(p1++) = *(p2++) + *(p3++);
    m.temporary_ = true;
    return m;
  }
}

template<typename DataType>
Matrix<DataType> Matrix<DataType>::operator+(
    const Matrix<DataType>& rhs) const {
  assert(num_col_ == rhs.num_col() && num_row_ == rhs.num_row());
  if (temporary_) {
    (const_cast<Matrix<DataType>* >(this))->operator+=(rhs);
    return *this;
  } else if (rhs.temporary_) {
    (const_cast<Matrix<DataType>& >(rhs)).operator+=(*this);
    return rhs;
  } else {
    Matrix<DataType> m(num_row_, num_col_);
    size_t endi = num_row_*num_col_;
    DataType *p1, *p2, *p3;
    p1 = m.raw_ptr();
    p2 = raw_ptr();
    p3 = rhs.raw_ptr();
    for (size_t i = 0; i < endi; i++)
      *(p1++) = *(p2++) + *(p3++);
    m.temporary_ = true;
    return m;
  }
}

template<typename DataType>
template<typename T>
Matrix<DataType>& Matrix<DataType>::operator+=(const Matrix<T>& rhs) {
  assert(num_col_ == rhs.num_col() && num_row_ == rhs.num_row());
  size_t endi = num_row_*num_col_;
  DataType *p2, *p3;
  p2 = raw_ptr();
  p3 = rhs.raw_ptr();
  for (size_t i = 0; i < endi; i++)
    *(p2++) += *(p3++);
  return *this;
}

template<typename DataType>
template<typename T>
Matrix<DataType> Matrix<DataType>::operator-(const Matrix<T>& rhs) const {
  assert(num_col_ == rhs.num_col() && num_row_ == rhs.num_row());
  if (temporary_) {
    (const_cast<Matrix<DataType>* >(this))->operator-=(rhs);
    return *this;
  } else {
    Matrix<DataType> m(num_row_, num_col_);
    size_t endi = num_row_*num_col_;
    DataType *p1, *p2, *p3;
    p1 = m.raw_ptr();
    p2 = raw_ptr();
    p3 = rhs.raw_ptr();
    for (size_t i = 0; i < endi; i++)
      *(p1++) = *(p2++) - *(p3++);
    temporary_ = true;
    return m;
  }
}

template<typename DataType>
Matrix<DataType> Matrix<DataType>::operator-(
    const Matrix<DataType>& rhs) const {
  assert(num_col_ == rhs.num_col() && num_row_ == rhs.num_row());
  if (temporary_) {
    (const_cast<Matrix<DataType>* >(this))->operator-=(rhs);
    return *this;
  } else if (rhs.temporary_) {
    size_t endi = num_row_*num_col_;
    DataType *p2, *p3;
    p2 = raw_ptr();
    p3 = rhs.raw_ptr();
    for (size_t i = 0; i < endi; i++) {
      *(p3) = *(p2++) - *(p3);
      p3++;
    }
    return rhs;
  } else {
    Matrix<DataType> m(num_row_, num_col_);
    size_t endi = num_row_*num_col_;
    DataType *p1, *p2, *p3;
    p1 = m.raw_ptr();
    p2 = raw_ptr();
    p3 = rhs.raw_ptr();
    for (size_t i = 0; i < endi; i++)
      *(p1++) = *(p2++) - *(p3++);
    m.temporary_ = true;
    return m;
  }
}

template<typename DataType>
template<typename T>
Matrix<DataType>& Matrix<DataType>::operator-=(const Matrix<T>& rhs) {
  assert(num_col_ == rhs.num_col() && num_row_ == rhs.num_row());
  size_t endi = num_row_*num_col_;
  DataType *p2, *p3;
  p2 = raw_ptr();
  p3 = rhs.raw_ptr();
  for (size_t i = 0; i < endi; i++)
    *(p2++) -= *(p3++);
  return *this;
}

template<typename DataType>
Matrix<DataType> Matrix<DataType>::operator*(double rhs) const {
  if (temporary_) {
    (const_cast<Matrix<DataType>* >(this))->operator*=(rhs);
    return *this;
  } else {
    Matrix<DataType> m(num_row_, num_col_);
    size_t endi = num_row_*num_col_;
    DataType *p1, *p2;
    p1 = m.raw_ptr();
    p2 = raw_ptr();
    for (size_t i = 0; i < endi; i++)
      *(p1++) = *(p2++)*rhs;
    m.temporary_ = true;
    return m;
  }
}

template<typename DataType>
Matrix<DataType> Matrix<DataType>::operator*(const CD &rhs) const {
  if (temporary_) {
    (const_cast<Matrix<DataType>* >(this))->operator*=(rhs);
    return *this;
  } else {
    Matrix<DataType> m(num_row_, num_col_);
    size_t endi = num_row_*num_col_;
    DataType *p1, *p2;
    p1 = m.raw_ptr();
    p2 = raw_ptr();
    for (size_t i = 0; i < endi; i++)
      *(p1++) = *(p2++)*rhs;
    m.temporary_ = true;
    return m;
  }
}

template<typename DataType>
Matrix<DataType> & Matrix<DataType>::operator*=(double rhs) {
  size_t endi = num_row_*num_col_;
  DataType *p2;
  p2 = raw_ptr();
  for (size_t i = 0; i < endi; i++)
    *(p2++) *= rhs;
  return *this;
}

template<typename DataType>
Matrix<DataType> & Matrix<DataType>::operator*=(const CD &rhs) {
  size_t endi = num_row_*num_col_;
  DataType *p2;
  p2 = raw_ptr();
  for (size_t i = 0; i < endi; i++)
    *(p2++) *= rhs;
  return *this;
}

template<typename DataType>
Matrix<DataType> Matrix<DataType>::operator*(
    const Matrix<DataType>& rhs) const {
  assert(num_col_ == rhs.num_row());
  size_t m, n, k1, k2, lda, ldb;
  m = num_row_;
  k1 = num_col_;
  k2 = rhs.num_row();
  n = rhs.num_col();
  lda = num_row_;
  ldb = rhs.num_row();

  Matrix<DataType> ma;

  ma.Resize(m, n);

  gemm<DataType>(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k1,
    raw_ptr(), lda, rhs.raw_ptr(), ldb, ma.raw_ptr(), m);

  return ma;
}

template<typename DataType>
Matrix<DataType>& Matrix<DataType>::operator*=(const Matrix<DataType>& rhs) {
  (*this) = (*this)*rhs;
  return (*this);
}

template<typename DataType>
template<typename T>
Matrix<DataType> Matrix<DataType>::operator/(const T &rhs) const {
  return operator*(1.0/rhs);
}

template<typename DataType>
template<typename T>
Matrix<DataType>& Matrix<DataType>::operator/=(const T &rhs) {
  return operator*=(1.0/rhs);
}

/* Matrix specific operations */
template<typename DataType>
Matrix<DataType> Matrix<DataType>::Transpose() const {
  if (temporary_ && num_row_ == num_col_) {
    size_t n;
    DataType temp;
    DataType *p0, *p1, *p2;
    p0 = raw_ptr();
    for (size_t j = 0; j < num_col_-1; j++) {
      n = j+1;
      p1 = p0 + j*num_row_ + n;
      p2 = p0 + j + num_col_*n;
      for (size_t i = n; i < num_row_; i++) {
        temp = *p1;
        *p1 = *p2;
        *p2 = temp;
        p1++;
        p2 += num_col_;
      }
    }
    const_cast<Matrix<DataType>*>(this)->Resize(num_col_, num_row_);
    return *this;
  } else if (temporary_ && (num_row_== 1 || num_col_ == 1)) {
    const_cast<Matrix<DataType>*>(this)->Resize(num_col_, num_row_);
    return *this;
  } else {
  Matrix<DataType> m(num_col_, num_row_);
  DataType *p0, *p1, *p2;
  p0 = m.raw_ptr();
  p2 = raw_ptr();
  for (size_t j = 0; j < num_col_; j++) {
    p1 = p0 + j;
    for (size_t i = 0; i < num_row_; i++) {
      *p1 = *(p2++);
      p1 += num_col_;
    }
  }
  m.set_temporary(true);
  return m;
  }
}

template<typename DataType>
Matrix<DataType> Matrix<DataType>::Adjoint() const {
  if (temporary_ && (num_row_ == num_col_ || num_row_ == 1 || num_col_ == 1)) {
    Conjugate();
    Transpose();
    return *this;
  } else {
    Matrix<DataType> m(num_col_, num_row_);
    DataType *p0, *p1, *p2;
    p0 = m.raw_ptr();
    p2 = raw_ptr();
    for (size_t j = 0; j < num_col_; j++) {
      p1 = p0 + j;
      for (size_t i = 0; i < num_row_; i++) {
        *p1 = mlcpp::Conjugate(*(p2++));
        p1 += num_col_;
      }
    }
    m.temporary_ = true;
    return m;
  }
}

template<typename DataType>
Matrix<DataType> Matrix<DataType>::Conjugate() const {
  if (temporary_) {
    size_t endi = num_row_*num_col_;
    DataType *p2;
    p2 = raw_ptr();
    for (size_t i = 0; i < endi; i++) {
      *(p2)=mlcpp::Conjugate(*(p2));
      p2++;
    }
    return *this;
  } else {
    Matrix<DataType> m(num_row_, num_col_);
    size_t endi = num_row_*num_col_;
    DataType *p1, *p2;
    p1 = m.raw_ptr();
    p2 = raw_ptr();
    for (size_t i = 0; i < endi; i++)
        *(p1++)=mlcpp::Conjugate(*(p2++));
    m.temporary_ = true;
    return m;
  }
}

template<typename DataType>
Matrix<DataType> Matrix<DataType>::Block(size_t r1, size_t r2,
                                         size_t c1, size_t c2) const {
  assert(r1 <= r2 && c1 <= c2);
  Matrix<DataType> m(r2-r1, c2-c1);
  DataType *p1, *p2;
  size_t n1, n2;
  n1 = num_row_;
  n2 = n1-r2+r1;
  p2 = m.raw_ptr();
  p1 = data_ptr_->data_ptr + r1+n1*c1;

  for (size_t j = c1; j < c2; j++) {
    for (size_t i = r1; i < r2; i++)
      *(p2++) = *(p1++);
    p1 += n2;
  }
  m.temporary_ = true;
  return m;
}

template<typename DataType>
Matrix<DataType>& Matrix<DataType>::Replace(size_t r, size_t c,
                                            const Matrix<DataType> &mt) {
  DataType *p1, *p2;
  size_t n1, n2, n3, n4;
  n1 = num_row_;
  p2 = mt.raw_ptr();
  n2 = mt.num_col();
  n3 = mt.num_row();
  n4 = n1-n3;
  p1 = data_ptr_->data_ptr + r+n1*c;

  for (size_t j = 0; j < n2; j++) {
    for (size_t i = 0; i < n3; i++)
      *(p1++) = *(p2++);
    p1 += n4;
  }

  return *this;
}

template<typename DataType>
inline Matrix<DataType> Matrix<DataType>::GetCol(size_t c) const {
  return Matrix<DataType>(raw_ptr_+c*num_row_, num_row_, 1);
}

template<typename DataType>
Matrix<DataType> Matrix<DataType>::GetRow(size_t r) const {
  Matrix<DataType> m(1, num_col_);
  DataType *p1, *p2;
  p2 = m.raw_ptr();
  p1 = raw_ptr_ + r;

  for (size_t j = 0; j < num_col_; j++) {
    *(p2++) = *(p1);
    p1 += num_row_;
  }
  m.temporary_ = true;
  return m;
}

template<typename DataType>
Matrix<DataType> & Matrix<DataType>::SwapCol(size_t c1, size_t c2) {
  DataType temp;
  DataType *p1, *p2;
  p1 = data_ptr_->data_ptr + c1*num_row_;
  p2 = data_ptr_->data_ptr + c2*num_row_;
  for (size_t i = 0; i < num_row_; i++) {
    temp = *p1;
    *(p1++) = *p2;
    *(p2++) = temp;
  }
  return *this;
}

template<typename DataType>
Matrix<DataType> & Matrix<DataType>::SwapRow(size_t r1, size_t r2) {
  DataType temp;
  DataType *p1, *p2;
  p1 = data_ptr_->data_ptr + r1;
  p2 = data_ptr_->data_ptr + r2;
  for (size_t i = 0; i < num_col_; i++) {
    temp = *p1;
    *p1 = *p2;
    *p2 = temp;
    p1 += num_row_;
    p2 += num_row_;
  }

  return *this;
}

/* solvers */
template<typename DataType>
template<typename T>
void Matrix<DataType>::EigenSolver(
    Matrix<Complex<T> > &e,
    Matrix<DataType> &vl,
    Matrix<DataType> &vr) const {
  assert(num_col_ == num_row_);
  Matrix<DataType> m(*this);

  vl.Resize(num_col_, num_col_);
  vr.Resize(num_col_, num_col_);
  e.Resize(num_col_, 1);

  char nl = &vl? 'V': 'N';
  char nr = &vr? 'V': 'N';

  geev<DataType>(nl, nr, num_col_, m.raw_ptr(), num_col_,
                 e.raw_ptr(), vl.raw_ptr(), num_col_,
                 vr.raw_ptr(), num_col_);
}

template<typename DataType>
template<typename T>
void Matrix<DataType>::RightEigenSolver(
    Matrix<Complex<T> > &e,
    Matrix<DataType> &vr) const {
  assert(num_col_ == num_row_);
  Matrix<DataType> m(*this);

  vr.Resize(num_col_, num_col_);
  e.Resize(num_col_, 1);

  char nl = 'N';
  char nr = 'V';

  geev<DataType>(nl, nr, num_col_, m.raw_ptr(), num_col_,
                 e.raw_ptr(), NULL, num_col_, vr.raw_ptr(), num_col_);
}

template<typename DataType>
template<typename T>
void Matrix<DataType>::LeftEigenSolver(
    Matrix<Complex<T> > &e,
    Matrix<DataType> &vl) const {
  assert(num_col_ == num_row_);
  Matrix<DataType> m(*this);

  vl.Resize(num_col_, num_col_);
  e.Resize(num_col_, 1);

  char nl = 'V';
  char nr = 'N';

  geev<DataType>(nl, nr, num_col_, m.raw_ptr(), num_col_,
                 e.raw_ptr(), vl.raw_ptr(), num_col_, NULL, num_col_);
}
/////////////////////////////////////////////////////////////////////////////

/**
 * Return the real part of a Matrix
 */
Matrix<float> Real(const Matrix<float> &ma) {
  Matrix<float> m(ma.num_row(), ma.num_col());
  size_t endi = ma.num_row()*ma.num_col();
  float *p1;
  float *p2;
  p1 = m.raw_ptr();
  p2 = ma.raw_ptr();
  for (size_t i = 0; i < endi; i++)
    *(p1++)=*(p2++);
  m.set_temporary(true);
  return m;
}

/**
 * Return the imaginary part of a Matrix
 */
Matrix<float> Imag(const Matrix<float> &ma) {
  Matrix<float> m(ma.num_row(), ma.num_col());
  m.Clear();
  m.set_temporary(true);
  return m;
}

/**
 * Return the real part of a Matrix
 */
Matrix<double> Real(const Matrix<double> &ma) {
  Matrix<double> m(ma.num_row(), ma.num_col());
  size_t endi = ma.num_row()*ma.num_col();
  double *p1;
  double *p2;
  p1 = m.raw_ptr();
  p2 = ma.raw_ptr();
  for (size_t i = 0; i < endi; i++)
    *(p1++)=*(p2++);
  m.set_temporary(true);
  return m;
}

/**
 * Return the imaginary part of a Matrix
 */
Matrix<double> Imag(const Matrix<double> &ma) {
  Matrix<double> m(ma.num_row(), ma.num_col());
  m.Clear();
  m.set_temporary(true);
  return m;
}

/**
 * Return the real part of a Matrix
 */
Matrix<float> Real(const Matrix<CS> &ma) {
  Matrix<float> m(ma.num_row(), ma.num_col());
  size_t endi = ma.num_row()*ma.num_col();
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
  size_t endi = ma.num_row()*ma.num_col();
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
  size_t endi = ma.num_row()*ma.num_col();
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
  size_t endi = ma.num_row()*ma.num_col();
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
template<typename DataType>
Matrix<DataType> Transpose(const Matrix<DataType> &m) {
  return m.Transpose();
}

/**
 * Return hermitian of a Matrix
 */
template<typename DataType>
Matrix<DataType> Adjoint(const Matrix<DataType> &m) {
  return m.Adjoint();
}

/**
 * Return conjugate of a Matrix
 */
template<typename DataType>
Matrix<DataType> Conjugate(const Matrix<DataType> &m) {
  return m.Conjugate();
}

/**
 * Multiply a double number and a Matrix
 */
template<typename DataType>
inline Matrix<DataType> operator*(double lhs,
                                 const Matrix<DataType> &ma) {
  return ma*lhs;
}

/**
 * Multiply a complex number and a Matrix
 */
template<typename DataType>
inline Matrix<DataType> operator*(const CD &lhs,
                                 const Matrix<DataType> &ma) {
  return ma*lhs;
}

/* stream operator overload */
template<typename DataType>
std::ostream& operator<< (std::ostream& out,
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
std::string OutputToString(const Matrix<DataType> &m) {
  std::ostringstream out;
  out << m;
  return out.str();
}

/////////////////////////////////////////////////////////////////////////////



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
  explicit IdentityMatrix(size_t n);
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
#endif  // MLCPP_MATRIX_H_  // NOLINT
