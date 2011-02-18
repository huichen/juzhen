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

#ifndef MLCPP_MATRIX_H_
#define MLCPP_MATRIX_H_
#include <mlcpp_complex.h>
#include <mlcpp_dataarray.h>
#include <sstream>

#include <mlcpp_adaptor.h>

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
  Matrix(const Matrix<T> &m);
 
  /**
   * Construct from another Matrix. All numbers are copied. 
   */
  Matrix(const Matrix<DataType> &m); 

  /**
   * Copy from another Matrix. New memory is allocated if it's necessary. 
   */
  Matrix<DataType>& operator= (const Matrix<DataType> &rhs);
 
  /**
   * Assign all numbers in the Matrix to be rhs.
   */
  Matrix<DataType>& operator= (const DataType rhs);
 
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
  size_t nrow() const;

  /**
   * Return number of columns that the Matrix has.
   */
  size_t ncol() const;

  /**
   * Returns the pointer of raw array. This is pretty useful when calling 
   * low level blas/lapack functions.
   */
  inline DataType * dataptr() const; 

  /**
   * Get a Matrix's temporary flag.
   */ 
  inline bool gettemporary() const;

  /**
   * Set a Matrix's temporary flag.
   */ 
  inline void settemporary(bool t);

  /**
   * Resizing or reshaping the Matrix. If the resized size is larger than
   * the original, new memory will be allocated and data will be copied to
   * the new position. 
   */
  void resize(size_t nr, size_t nc); 

  /**
   * Set all elements to be zero.
   */
  inline void clear(); 

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
  const Matrix<DataType> operator+(const Matrix<T>& rhs) const;
 
  /** 
   * Add two matrices
   */
  const Matrix<DataType> operator+(const Matrix<DataType>& rhs) const;
 
  /** 
   * Add one matrices upto another Matrix. 
   */
  template<typename T> 
  Matrix<DataType>& operator+=(const Matrix<T>& rhs);
 
  /**
   * Subtract one Matrix and another Matrix.
   */
  template<typename T> 
  const Matrix<DataType> operator-(const Matrix<T>& rhs) const;

  /**
   * Subtract one Matrix and another Matrix.
   */
  const Matrix<DataType> operator-(const Matrix<DataType>& rhs) const;

  /**
   * Subtract one Matrix from another Matrix.
   */
  template<typename T> 
  Matrix<DataType>& operator-=(const Matrix<T>& rhs);

  /**
   * Multiply each element in the Matrix by a real number and
   * get a new Matrix.
   */
  const Matrix<DataType> operator*(const double rhs) const;

  /**
   * Multiply each element in the Matrix by a complex number and 
   * get a new Matrix.
   */
  const Matrix<DataType> operator*(const CD rhs) const;

  /**
   * Multiply each element in the Matrix by a real number.
   */
  Matrix<DataType> & operator*=(const double rhs);

  /**
   * Multiply each element in the Matrix by a complex number.
   */
  Matrix<DataType> & operator*=(const CD rhs);
 
  /**
   * Multiply two matrices and get a new Matrix. Number of columns of the first
   * Matrix must equals to the number of columns of the second Matrix.
   */
  const Matrix<DataType> operator*(const Matrix<DataType>& rhs) const;

  /**
   * Multiply the Matrix with a new Matrix. Number of columns of the first
   * Matrix must equals to the number of columns of the second Matrix.
   */
  Matrix<DataType>& operator*=(const Matrix<DataType>& rhs);
 
  /**
   * Divide each element in the Matrix by a constant to get a new Matrix.
   */
  template<typename T> 
  const Matrix<DataType> operator/(const T rhs) const;

  /**
   * Divide each element in the Matrix by a constant.
   */
  template<typename T> 
  Matrix<DataType>& operator/=(const T rhs);

  /*************************************
   *  Matrix specific operations 
   ************************************/
  /**
   * Get the Matrix's transposition. 
   */
  Matrix<DataType> trans()const;

  /**
   * Get the Matrix's hermitian.  
   */
  Matrix<DataType> herm() const;

  /**
   * Get the Matrix's conjugate.  
   */
  Matrix<DataType> conj()const;

  /**
   * Get a sub-block of the Matrix between r1th (containing) row and r2th 
   * (not containing) row, and between c1th (containing) column and c2th 
   * (not containing) column.  
   */
  Matrix<DataType> block(size_t r1, size_t r2, size_t c1, size_t c2) const;

  /**
   * Replace a sub-block of the Matrix with another Matrix value. The
   * sub-block starts at (r-th row, c-th column) at its upper-left corner,
   * and is replace by values from mt's (0,0) position. 
   */
  const Matrix<DataType> &replace(size_t r, size_t c, const Matrix<DataType> &mt) const;

  /**
   * Return the Matrix's c-th column.
   */
  inline Matrix<DataType> col(size_t c) const;

  /**
   * Return the Matrix's r-th row.
   */
  Matrix<DataType> row(size_t r) const;

  /**
   * Swap two cols.
   */
  Matrix<DataType> &swapcol(size_t c1, size_t c2);

  /**
   * Swap two rows.
   */
  Matrix<DataType> &swaprow(size_t r1, size_t r2);

  /**
   * Solve eigen system and put eigen values into e, corresponding 
   * left-eigen vectors into vl and corresponding right-eigen vectors
   * into vr. vl and vr will be square matrices.
   * 
   * If you are linking to basic blas library only, this function is not
   * implemented yet. 
   */
  template<typename T>
  void eigen(Matrix<Complex<T> > &e, Matrix<DataType> &vl, Matrix<DataType> &vr) const;

  /**
   * Solve eigen system and put eigen values into e, corresponding 
   * left-eigen vectors into vl. vl will be square Matrix.
   *
   * If you are linking to basic blas library only, this function is not
   * implemented yet. 
   */
  template<typename T>
  void reigen(Matrix<Complex<T> > &e, Matrix<DataType> &vr) const;

  /**
   * Solve eigen system and put eigen values into e, corresponding 
   * right-eigen vectors into vr. vr will be square Matrix.
   *
   * If you are linking to basic blas library only, this function is not
   * implemented yet. 
   */
  template<typename T>
  void leigen(Matrix<Complex<T> > &e, Matrix<DataType> &vl) const;

  typedef std::auto_ptr<DataArray<DataType> > DataPtr;

protected:
  /**
   * Pointer to raw array. 
   */
  DataType * m_rawptr; 

  /**
   * Auto pointer to m_rawptr. 
   */
  DataPtr m_data; 

  /**
   * Number of columns.
   */
  size_t m_ncol;

  /**
   * Number of rows.
   */
  size_t m_nrow;

  /**
   * If the Matrix is temporary. 
   */
  bool m_temporary;

};

typedef Matrix<float> smatrix;
typedef Matrix<double> dmatrix;
typedef Matrix<CS> cmatrix;
typedef Matrix<CD> zmatrix;

template<typename DataType> 
Matrix<DataType>::Matrix() : m_ncol(0), m_nrow(0) {
  m_temporary = false;
  m_rawptr = NULL;
} 

template<typename DataType> 
Matrix<DataType>::Matrix(size_t nr, size_t nc) : m_ncol(nc), m_nrow(nr) {
  m_data = DataPtr(new DataArray<DataType>(m_ncol*m_nrow));
  m_temporary = false;
  m_rawptr = m_data->m_data;
} 

template<typename DataType> 
template<typename T> 
Matrix<DataType>::Matrix(const T *data, size_t nr, size_t nc) 
  : m_ncol(nc), m_nrow(nr) {
  m_data =DataPtr(new DataArray<DataType>(data, nr*nc));
  m_temporary = false;
  m_rawptr = m_data->m_data;
} 

template<typename DataType> 
template<typename T> 
Matrix<DataType>::Matrix(const Matrix<T> &m) 
  : m_ncol(m.ncol()), m_nrow(m.nrow()) { 
  m_data =DataPtr(
    new DataArray<DataType>(m.dataptr(), m.nrow()*m.ncol())
  ); 
  m_temporary = false;
  m_rawptr = m_data->m_data;
} 

template<typename DataType> 
Matrix<DataType>::Matrix(const Matrix<DataType> &m) 
  : m_ncol(m.ncol()), m_nrow(m.nrow()) { 
  if (m.m_temporary) {
    m_data = (const_cast<Matrix<DataType>& >(m)).m_data;
    m_temporary = true;
  } else {
    m_data =DataPtr(new DataArray<DataType>(*(m.m_data)));
    m_temporary = false;
  }
  m_rawptr = m_data->m_data;
}

template<typename DataType> 
Matrix<DataType>& 
Matrix<DataType>::operator= (const Matrix<DataType> &rhs) {
  if (&rhs==this) return *this;
  if (rhs.m_temporary) {
    m_data = (const_cast<Matrix<DataType>& >(rhs)).m_data;
    m_ncol = rhs.ncol();
    m_nrow = rhs.nrow();
    m_rawptr = m_data->m_data;
    return *this;
  }
  if ((m_rawptr == NULL) || m_data->m_size < rhs.ncol()*rhs.nrow()) 
    m_data =DataPtr(new DataArray<DataType>(*(rhs.m_data)));
  else 
    memcpy(m_data->m_data, rhs.dataptr(), rhs.ncol()*rhs.nrow()*sizeof(DataType));
  m_ncol = rhs.ncol();
  m_nrow = rhs.nrow();
  m_rawptr = m_data->m_data;
  return *this;
} 

template<typename DataType> 
Matrix<DataType>& 
Matrix<DataType>::operator= (const DataType rhs) {
  size_t endi = m_ncol*m_nrow; 
  DataType *p = m_data->m_data;
  for (size_t i=0; i<endi; i++) *(p++)=rhs;
  return *this;
} 

template<typename DataType> 
bool Matrix<DataType>::operator== (const Matrix<DataType> &rhs) {
  if (&rhs==this) return true;
  if (m_nrow != rhs.nrow() || m_ncol != rhs.ncol()) return false;
  size_t endi = m_ncol*m_nrow; 
  DataType *p1 = m_data->m_data;
  DataType *p2 = rhs.dataptr();
  for (size_t i=0; i<endi; i++)
    if (*(p1++)!=*(p2++)) return false;
  return true;
} 

template<typename DataType> 
bool Matrix<DataType>::operator!= (const Matrix<DataType> &rhs) {
  return !(operator==(rhs));
}

/* interface to private data */
template<typename DataType> 
size_t Matrix<DataType>::nrow() const {
  return m_nrow;
}

template<typename DataType> 
size_t Matrix<DataType>::ncol() const {
  return m_ncol;
}

template<typename DataType> 
inline DataType * Matrix<DataType>::dataptr() const {
//  return m_data->m_data;
  return m_rawptr;
}

/* resizing/reshaping Matrix */
template<typename DataType> 
void Matrix<DataType>::resize(size_t nr, size_t nc) {
  if (!m_rawptr || m_data->m_size < nc * nr) {
    DataPtr 
      newData(new DataArray<DataType>(*m_data, nr*nc));
    m_data = newData;
    m_rawptr = m_data->m_data;
  }
  m_ncol = nc;
  m_nrow = nr;
} 

template<typename DataType>
inline bool Matrix<DataType>::gettemporary() const {
  return m_temporary;
}

template<typename DataType>
inline void Matrix<DataType>::settemporary(bool t) {
  m_temporary = t;
};

template<typename DataType> 
inline void Matrix<DataType>::clear() {
  if (m_data->m_data)
    memset(m_data->m_data, 0, m_data->m_size*sizeof(DataType));
}

/* operator overloading */
template<typename DataType> 
inline DataType& Matrix<DataType>::operator()(size_t i, size_t j) {
  assert(i<m_nrow && j<m_ncol);
  return m_rawptr[j*m_nrow + i];
}

template<typename DataType> 
inline DataType Matrix<DataType>::operator()(size_t i, size_t j) const {
  assert(i<m_nrow && j<m_ncol);
  return m_rawptr[j*m_nrow + i];
}

template<typename DataType> 
inline DataType& Matrix<DataType>::operator()(size_t i) {
  return m_rawptr[i];
}

template<typename DataType> 
inline DataType& Matrix<DataType>::operator()(size_t i) const {
  return m_rawptr[i];
}

template<typename DataType> 
inline DataType& Matrix<DataType>::operator[](size_t i) {
  return m_rawptr[i];
}

template<typename DataType> 
inline DataType& Matrix<DataType>::operator[](size_t i) const {
  return m_rawptr[i];
}

// arithmetic

template<typename DataType> 
template<typename T> 
const Matrix<DataType> 
Matrix<DataType>::operator+(const Matrix<T>& rhs) const {
  assert(m_ncol == rhs.ncol() && m_nrow == rhs.nrow());
  if (m_temporary) {
    (const_cast<Matrix<DataType>* >(this))->operator+=(rhs);
    return *this;
  } else {
    Matrix<DataType> m(m_nrow, m_ncol);
    size_t endi = m_nrow*m_ncol;
    DataType *p1, *p2, *p3;
    p1=m.dataptr();
    p2=dataptr();
    p3=rhs.dataptr();
    for (size_t i=0; i<endi; i++)
      *(p1++) = *(p2++) + *(p3++);
    m.m_temporary = true;
    return m;
  }
}

template<typename DataType> 
const Matrix<DataType> 
Matrix<DataType>::operator+(const Matrix<DataType>& rhs) const {
  assert(m_ncol == rhs.ncol() && m_nrow == rhs.nrow());
  if (m_temporary) {
    (const_cast<Matrix<DataType>* >(this))->operator+=(rhs);
    return *this;
  } else if (rhs.m_temporary) {
    (const_cast<Matrix<DataType>& >(rhs)).operator+=(*this);
    return rhs;
  } else {
    Matrix<DataType> m(m_nrow, m_ncol);
    size_t endi = m_nrow*m_ncol;
    DataType *p1, *p2, *p3;
    p1=m.dataptr();
    p2=dataptr();
    p3=rhs.dataptr();
    for (size_t i=0; i<endi; i++)
      *(p1++) = *(p2++) + *(p3++);
    m.m_temporary = true;
    return m;
  }
}

template<typename DataType> 
template<typename T> 
Matrix<DataType>& Matrix<DataType>::operator+=(const Matrix<T>& rhs) {
  assert(m_ncol == rhs.ncol() && m_nrow == rhs.nrow());
  size_t endi = m_nrow*m_ncol;
  DataType *p2, *p3;
  p2=dataptr();
  p3=rhs.dataptr();
  for (size_t i=0; i<endi; i++)
    *(p2++) += *(p3++);
  return *this;
}

template<typename DataType> 
template<typename T>
const Matrix<DataType> 
Matrix<DataType>::operator-(const Matrix<T>& rhs) const {
  assert(m_ncol == rhs.ncol() && m_nrow == rhs.nrow());
  if (m_temporary) {
    (const_cast<Matrix<DataType>* >(this))->operator-=(rhs);
    return *this;
  } else {
    Matrix<DataType> m(m_nrow, m_ncol);
    size_t endi = m_nrow*m_ncol;
    DataType *p1, *p2, *p3;
    p1=m.dataptr();
    p2=dataptr();
    p3=rhs.dataptr();
    for (size_t i=0; i<endi; i++)
      *(p1++) = *(p2++) - *(p3++);
    m_temporary = true;
    return m;
  }
}

template<typename DataType> 
const Matrix<DataType> 
Matrix<DataType>::operator-(const Matrix<DataType>& rhs) const {
  assert(m_ncol == rhs.ncol() && m_nrow == rhs.nrow());
  if (m_temporary) {
    (const_cast<Matrix<DataType>* >(this))->operator-=(rhs);
    return *this;
  } else if (rhs.m_temporary) {
    size_t endi = m_nrow*m_ncol;
    DataType *p2, *p3;
    p2=dataptr();
    p3=rhs.dataptr();
    for (size_t i=0; i<endi; i++) {
      *(p3) = *(p2++) - *(p3);
      p3++;
    }
    return rhs;
  } else {
    Matrix<DataType> m(m_nrow, m_ncol);
    size_t endi = m_nrow*m_ncol;
    DataType *p1, *p2, *p3;
    p1=m.dataptr();
    p2=dataptr();
    p3=rhs.dataptr();
    for (size_t i=0; i<endi; i++)
      *(p1++) = *(p2++) - *(p3++);
    m.m_temporary = true;
    return m;
  }
}

template<typename DataType> 
template<typename T> 
Matrix<DataType>& Matrix<DataType>::operator-=(const Matrix<T>& rhs) {
  assert(m_ncol == rhs.ncol() && m_nrow == rhs.nrow());
  size_t endi = m_nrow*m_ncol;
  DataType *p2, *p3;
  p2=dataptr();
  p3=rhs.dataptr();
  for (size_t i=0; i<endi; i++)
    *(p2++) -= *(p3++);
  return *this;
}

template<typename DataType> 
const Matrix<DataType> Matrix<DataType>::operator*(const double rhs) const {
  if (m_temporary) {
    (const_cast<Matrix<DataType>* >(this))->operator*=(rhs);
    return *this;
  } else {
    Matrix<DataType> m(m_nrow, m_ncol);
    size_t endi = m_nrow*m_ncol;
    DataType *p1, *p2;
    p1=m.dataptr();
    p2=dataptr();
    for (size_t i=0; i<endi; i++)
      *(p1++) = *(p2++)*rhs; 
    m.m_temporary = true;
    return m;
  }
}

template<typename DataType> 
const Matrix<DataType> Matrix<DataType>::operator*(const CD rhs) const {
  if (m_temporary) {
    (const_cast<Matrix<DataType>* >(this))->operator*=(rhs);
    return *this;
  } else {
    Matrix<DataType> m(m_nrow, m_ncol);
    size_t endi = m_nrow*m_ncol;
    DataType *p1, *p2;
    p1=m.dataptr();
    p2=dataptr();
    for (size_t i=0; i<endi; i++)
      *(p1++) = *(p2++)*rhs; 
    m.m_temporary = true;
    return m;
  }
}

template<typename DataType> 
Matrix<DataType> & Matrix<DataType>::operator*=(const double rhs) {
  size_t endi = m_nrow*m_ncol;
  DataType *p2;
  p2=dataptr();
  for (size_t i=0; i<endi; i++)
    *(p2++) *= rhs; 
  return *this;
}

template<typename DataType> 
Matrix<DataType> & Matrix<DataType>::operator*=(const CD rhs) {
  size_t endi = m_nrow*m_ncol;
  DataType *p2;
  p2=dataptr();
  for (size_t i=0; i<endi; i++)
    *(p2++) *= rhs; 
  return *this;
}
 
template<typename DataType> 
const Matrix<DataType> 
Matrix<DataType>::operator*(const Matrix<DataType>& rhs) const {
  assert (m_ncol == rhs.nrow());
  size_t m,n,k1,k2,lda,ldb;
  m = m_nrow;
  k1 = m_ncol;
  k2 = rhs.nrow();
  n = rhs.ncol();
  lda = m_nrow;
  ldb = rhs.nrow();

  Matrix<DataType> ma;

  ma.resize(m,n);

  gemm<DataType>(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k1, 
    dataptr(), lda, rhs.dataptr(), ldb, ma.dataptr(), m);

  return ma;
} 

template<typename DataType> 
Matrix<DataType>& Matrix<DataType>::operator*=(const Matrix<DataType>& rhs) {
  (*this) = (*this)*rhs; 
  return (*this);
}
 
template<typename DataType> 
template<typename T> 
const Matrix<DataType> Matrix<DataType>::operator/(const T rhs) const {
  return operator*(1.0/rhs);
}

template<typename DataType> 
template<typename T> 
Matrix<DataType>& Matrix<DataType>::operator/=(const T rhs) {
  return operator*=(1.0/rhs);
}

/* Matrix specific operations */
template<typename DataType> 
Matrix<DataType> Matrix<DataType>::trans() const {
  if (m_temporary && m_nrow == m_ncol) {
    size_t n;
    DataType temp;
    DataType *p0, *p1, *p2;
    p0=dataptr();
    for (size_t j=0; j<m_ncol-1; j++) { 
      n = j+1;
      p1 = p0 + j*m_nrow + n;
      p2 = p0 + j + m_ncol*n; 
      for (size_t i=n; i<m_nrow; i++) { 
        temp = *p1;
        *p1 = *p2; 
        *p2 = temp;
        p1++;
        p2+=m_ncol;
      }
    }
    const_cast<Matrix<DataType>*>(this)->resize(m_ncol, m_nrow);
    return *this;
  } else if (m_temporary && (m_nrow==1 || m_ncol ==1) ){
    const_cast<Matrix<DataType>*>(this)->resize(m_ncol, m_nrow);
    return *this;
  } else {
  Matrix<DataType> m(m_ncol, m_nrow); 
  DataType *p0, *p1, *p2;
  p0=m.dataptr();
  p2=dataptr();
  for (size_t j=0; j<m_ncol; j++) { 
    p1 = p0 + j;
    for (size_t i=0; i<m_nrow; i++) { 
      *p1=*(p2++);
      p1+=m_ncol;
    }
  }
  m.settemporary(true);
  return m;
  }
} 

template<typename DataType> 
Matrix<DataType> Matrix<DataType>::herm() const {
  if (m_temporary && (m_nrow == m_ncol || m_nrow == 1 || m_ncol == 1)) {
    conj();
    trans();
    return *this;
  } else {
    Matrix<DataType> m(m_ncol, m_nrow); 
    DataType *p0, *p1, *p2;
    p0=m.dataptr();
    p2=dataptr();
    for (size_t j=0; j<m_ncol; j++) { 
      p1 = p0 + j;
      for (size_t i=0; i<m_nrow; i++) { 
        *p1=mlcpp::conj(*(p2++));
        p1+=m_ncol;
      }
    }
    m.m_temporary = true;
    return m;
  }
} 

template<typename DataType> 
Matrix<DataType> Matrix<DataType>::conj() const {
  if (m_temporary) {
    size_t endi = m_nrow*m_ncol;
    DataType *p2;
    p2=dataptr();
    for (size_t i=0; i<endi; i++) { 
      *(p2)=mlcpp::conj(*(p2));
      p2++;
    }
    return *this;
  } else {
    Matrix<DataType> m(m_nrow, m_ncol); 
    size_t endi = m_nrow*m_ncol;
    DataType *p1, *p2;
    p1=m.dataptr();
    p2=dataptr();
    for (size_t i=0; i<endi; i++) 
        *(p1++)=mlcpp::conj(*(p2++));
    m.m_temporary = true;
    return m;
  }
} 

template<typename DataType> 
Matrix<DataType> Matrix<DataType>::block(size_t r1, size_t r2, 
                                         size_t c1, size_t c2) const {
  assert(r1<=r2 && c1<=c2);
  Matrix<DataType> m(r2-r1, c2-c1);
  DataType *p1, *p2;
  size_t n1, n2;
  n1 = m_nrow;
  n2 = n1-r2+r1;
  p2 = m.dataptr(); 
  p1 = m_data->m_data + r1+n1*c1;
  
  for (size_t j=c1; j<c2; j++) { 
    for (size_t i=r1; i<r2; i++) 
      *(p2++) = *(p1++);
    p1 += n2;
  }
  m.m_temporary = true;
  return m;
}

template<typename DataType> 
const Matrix<DataType> & Matrix<DataType>::replace(size_t r, size_t c, 
                                            const Matrix<DataType> &mt) const {
  DataType *p1, *p2;
  size_t n1, n2, n3, n4;
  n1 = m_nrow;
  p2 = mt.dataptr(); 
  n2 = mt.ncol();
  n3 = mt.nrow();
  n4 = n1-n3;
  p1 = m_data->m_data + r+n1*c;
  
  for (size_t j=0; j<n2; j++) { 
    for (size_t i=0; i<n3; i++) 
      *(p1++) = *(p2++);
    p1 += n4;
  }

  return *this;
}

template<typename DataType> 
inline Matrix<DataType> Matrix<DataType>::col(size_t c) const {
  return Matrix<DataType>(m_rawptr+c*m_nrow, m_nrow, 1);
}

template<typename DataType> 
Matrix<DataType> Matrix<DataType>::row(size_t r) const {
  Matrix<DataType> m(1, m_ncol);
  DataType *p1, *p2;
  p2 = m.dataptr(); 
  p1 = m_rawptr + r;
  
  for (size_t j=0; j<m_ncol; j++) { 
    *(p2++) = *(p1);
    p1 += m_nrow;
  }
  m.m_temporary = true;
  return m;
}

template<typename DataType> 
Matrix<DataType> & Matrix<DataType>::swapcol(size_t c1, size_t c2) {
  DataType temp;
  DataType *p1, *p2;
  p1 = m_data->m_data + c1*m_nrow;
  p2 = m_data->m_data + c2*m_nrow;
  for(size_t i=0; i<m_nrow; i++) {
    temp = *p1;
    *(p1++) = *p2;
    *(p2++) = temp;
  }
  return *this;
}

template<typename DataType> 
Matrix<DataType> & Matrix<DataType>::swaprow(size_t r1, size_t r2) {
  DataType temp;
  DataType *p1, *p2;
  p1 = m_data->m_data + r1;
  p2 = m_data->m_data + r2;
  for(size_t i=0; i<m_ncol; i++) {
    temp = *p1;
    *p1 = *p2;
    *p2 = temp;
    p1 += m_nrow;
    p2 += m_nrow;
  }

  return *this;
}

/* solvers */
template<typename DataType> 
template<typename T> 
void Matrix<DataType>::eigen(Matrix<Complex<T> > &e, Matrix<DataType> &vl, 
                             Matrix<DataType> &vr) const {
  assert (m_ncol == m_nrow);
  Matrix<DataType> m(*this);

  vl.resize(m_ncol,m_ncol);
  vr.resize(m_ncol,m_ncol);
  e.resize(m_ncol, 1);
  
  char nl = &vl? 'V': 'N';
  char nr = &vr? 'V': 'N';

  geev<DataType>(nl, nr, m_ncol, m.dataptr(), m_ncol, 
                 e.dataptr(), vl.dataptr(), m_ncol, 
                 vr.dataptr(), m_ncol);
} 

template<typename DataType> 
template<typename T> 
void Matrix<DataType>::reigen(Matrix<Complex<T> > &e, Matrix<DataType> &vr) const {
  assert (m_ncol == m_nrow);
  Matrix<DataType> m(*this);

  vr.resize(m_ncol,m_ncol);
  e.resize(m_ncol, 1);
  
  char nl = 'N';
  char nr = 'V';

  geev<DataType>(nl, nr, m_ncol, m.dataptr(), m_ncol, 
                 e.dataptr(), NULL, m_ncol, vr.dataptr(), m_ncol); 
} 

template<typename DataType> 
template<typename T> 
void Matrix<DataType>::leigen(Matrix<Complex<T> > &e, Matrix<DataType> &vl) const {
  assert (m_ncol == m_nrow);
  Matrix<DataType> m(*this);

  vl.resize(m_ncol,m_ncol);
  e.resize(m_ncol, 1);
  
  char nl = 'V';
  char nr = 'N';

  geev<DataType>(nl, nr, m_ncol, m.dataptr(), m_ncol, 
                 e.dataptr(), vl.dataptr(), m_ncol, NULL, m_ncol); 
} 
/////////////////////////////////////////////////////////////////////////////

/** 
 * Return the real part of a Matrix 
 */
Matrix<float> real(const Matrix<float> &ma) {
  Matrix<float> m(ma.nrow(), ma.ncol());
  size_t endi = ma.nrow()*ma.ncol();
  float *p1;
  float *p2;
  p1=m.dataptr();
  p2=ma.dataptr();
  for (size_t i=0; i<endi; i++) 
    *(p1++)=*(p2++);
  m.settemporary(true);
  return m;
} 

/** 
 * Return the imaginary part of a Matrix 
 */
Matrix<float> imag(const Matrix<float> &ma) {
  Matrix<float> m(ma.nrow(), ma.ncol());
  m.clear();
  m.settemporary(true);
  return m;
} 

/** 
 * Return the real part of a Matrix 
 */
Matrix<double> real(const Matrix<double> &ma) {
  Matrix<double> m(ma.nrow(), ma.ncol());
  size_t endi = ma.nrow()*ma.ncol();
  double *p1;
  double *p2;
  p1=m.dataptr();
  p2=ma.dataptr();
  for (size_t i=0; i<endi; i++) 
    *(p1++)=*(p2++);
  m.settemporary(true);
  return m;
} 

/** 
 * Return the imaginary part of a Matrix 
 */
Matrix<double> imag(const Matrix<double> &ma) {
  Matrix<double> m(ma.nrow(), ma.ncol());
  m.clear();
  m.settemporary(true);
  return m;
} 

/** 
 * Return the real part of a Matrix 
 */
Matrix<float> real(const Matrix<CS> &ma) {
  Matrix<float> m(ma.nrow(), ma.ncol());
  size_t endi = ma.nrow()*ma.ncol();
  float *p1;
  CS *p2;
  p1=m.dataptr();
  p2=ma.dataptr();
  for (size_t i=0; i<endi; i++) 
    *(p1++)=(p2++)->real;
  m.settemporary(true);
  return m;
} 

/** 
 * Return the imaginary part of a Matrix 
 */
Matrix<float> imag(const Matrix<CS> &ma) {
  Matrix<float> m(ma.nrow(), ma.ncol());
  size_t endi = ma.nrow()*ma.ncol();
  float *p1;
  CS *p2;
  p1=m.dataptr();
  p2=ma.dataptr();
  for (size_t i=0; i<endi; i++) 
    *(p1++)=(p2++)->imag;
  m.settemporary(true);
  return m;
} 

/** 
 * Return the real part of a Matrix 
 */
Matrix<double> real(const Matrix<CD> &ma) {
  Matrix<double> m(ma.nrow(), ma.ncol());
  size_t endi = ma.nrow()*ma.ncol();
  double *p1;
  CD *p2;
  p1=m.dataptr();
  p2=ma.dataptr();
  for (size_t i=0; i<endi; i++) 
    *(p1++)=(p2++)->real;
  m.settemporary(true);
  return m;
} 

/** 
 * Return the imaginary part of a Matrix 
 */
Matrix<double> imag(const Matrix<CD> &ma) {
  Matrix<double> m(ma.nrow(), ma.ncol());
  size_t endi = ma.nrow()*ma.ncol();
  double *p1;
  CD *p2;
  p1=m.dataptr();
  p2=ma.dataptr();
  for (size_t i=0; i<endi; i++) 
    *(p1++)=(p2++)->imag;
  m.settemporary(true);
  return m;
} 


/** 
 * Return transpose of a Matrix 
 */
template<typename DataType> 
Matrix<DataType> trans(const Matrix<DataType> &m) {
  return m.trans();
}

/** 
 * Return hermitian of a Matrix 
 */
template<typename DataType> 
Matrix<DataType> herm(const Matrix<DataType> &m) {
  return m.herm();
}

/** 
 * Return conjugate of a Matrix 
 */
template<typename DataType> 
Matrix<DataType> conj(const Matrix<DataType> &m) {
  return m.conj();
}

/** 
 * Multiply a double number and a Matrix
 */
template<typename DataType> 
inline const Matrix<DataType> operator*(const double lhs, 
                                 const Matrix<DataType> &ma) {
  return ma*lhs;
}

/** 
 * Multiply a complex number and a Matrix
 */
template<typename DataType> 
inline const Matrix<DataType> operator*(const CD lhs, 
                                 const Matrix<DataType> &ma) {
  return ma*lhs;
}

/* stream operator overload */
template<typename DataType> 
std::ostream& operator<< (std::ostream& out, 
                          const Matrix<DataType> &m) {
  for (size_t i=0; i<m.nrow(); i++) {
    for (size_t j=0; j<m.ncol(); j++)
      out << m(i,j) << " ";
    if (i!=m.nrow()-1) out << std::endl;
  }
  return out; 
}

/**
 * Get print form of a Matrix.
 */
template<typename DataType> 
std::string toString(const Matrix<DataType> &m) {
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
class IdMatrix : public Matrix<DataType> {
public:
  /**
   * Construct an identity Matrix of n x n. 
   */
  IdMatrix(size_t n);
};

typedef IdMatrix<float> identity_smatrix;
typedef IdMatrix<double> identity_dmatrix;
typedef IdMatrix<CS> identity_cmatrix;
typedef IdMatrix<CD> identity_zmatrix;

template<typename DataType> 
IdMatrix<DataType>::IdMatrix(size_t n) {
  Matrix<DataType>::m_data = typename Matrix<DataType>::DataPtr(new DataArray<DataType>(n*n));
  Matrix<DataType>::m_ncol = n;
  Matrix<DataType>::m_nrow = n;
  Matrix<DataType>::m_rawptr = Matrix<DataType>::m_data->m_data;
  
  size_t endi=n*n;
  DataType *p=Matrix<DataType>::dataptr();
  for (size_t i=0; i<endi; i++) 
    *(p++) = 0;

  p=Matrix<DataType>::dataptr();
  for (size_t i=0; i<n; i++) 
    *(p+i*(n+1)) = 1.;

}
/////////////////////////////////////////////////////////////////////////////


}
#endif
