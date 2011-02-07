#ifndef MLCPP_MATRIX_HPP
#define MLCPP_MATRIX_HPP
#include <mlcpp_complex.hpp>
#include <mlcpp_dataarray.hpp>
#include <sstream>

#include <mlcpp_wrapper.hpp>

namespace mlcpp {

/////////////////////////////////////////////////////////////////////////////
/** 
 * The basic matrix class that all variant matrix classes are derived from
 */
template<typename DataType> 
class matrix {
public:
  /**
   * Auto pointer to resource. 
   */
  typename DataPtr<DataType>::Type m_data; 

  /**
   * If the matrix is temporary. 
   */
  bool m_temporary;

  /**
   * Deconstructor. 
   */
  ~matrix() {}

  /**
   * Default constructor. 
   */
  matrix();

  /**
   * Construct an matrix of nr rows and nc columns. Memory is allocated but
   * matrix's items have undetermined values. 
   */
  matrix(size_t nr, size_t nc);
 
  /**
   * Construct an matrix of nr rows and nc columns from a raw array. nr*nc
   * numbers will be allocated in memory. If data's size is larger than nr*nc 
   * only first nr*nc numbers are used in column-major order. If data's 
   * size is smaller than nr*nc, all numbers from data will be used and the 
   * matrix's rests number have undetermined values. 
   */
  template<typename T> 
  matrix(const T *data, size_t nr, size_t nc); 

  /**
   * Construct from another matrix. All numbers are copied. 
   */
  template<typename T> 
  matrix(const matrix<T> &m);
 
  /**
   * Construct from another matrix. All numbers are copied. 
   */
  matrix(const matrix<DataType> &m); 

  /**
   * Copy from another matrix. New memory is allocated if it's necessary. 
   */
  matrix<DataType>& operator= (const matrix<DataType> &rhs);
 
  /**
   * Assign all numbers in the matrix to be rhs.
   */
  matrix<DataType>& operator= (const DataType rhs);
 
  /**
   * Check if two matrices are equal. The two matrices must have the same
   * numbers of rows and columns if it returns true. 
   */
  bool operator== (const matrix<DataType> &rhs); 

  /**
   * Check if two matrices are not equal. Returns false if the two matrices 
   * have different numbers of rows or columns. 
   */
  bool operator!= (const matrix<DataType> &rhs); 

  /***********************************************
   *  interface to private data 
   **********************************************/

  /**
   * Return number of rows that the matrix has.
   */
  size_t nrow() const;

  /**
   * Return number of columns that the matrix has.
   */
  size_t ncol() const;

  /**
   * Returns the pointer of raw array. This is pretty useful when calling 
   * low level blas/lapack functions.
   */
  DataType * getDataPtr() const; 

  /**
   * Return the reference to the matrix's DataArray.
   */
  DataArray<DataType>& getDataArray() const; 

  /**
   * Resizing or reshaping the matrix. If the resized size is larger than
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
  const matrix<DataType> operator+(const matrix<T>& rhs) const;
 
  /** 
   * Add one matrices upto another matrix. 
   */
  template<typename T> matrix<DataType>& operator+=(const matrix<T>& rhs);
 
  /**
   * Subtract one matrix and another matrix.
   */
  template<typename T> const matrix<DataType> operator-(const matrix<T>& rhs) const;

  /**
   * Subtract one matrix from another matrix.
   */
  template<typename T> matrix<DataType>& operator-=(const matrix<T>& rhs);

  /**
   * Multiply each element in the matrix by a real number and
   * get a new matrix.
   */
  const matrix<DataType> operator*(const double rhs) const;

  /**
   * Multiply each element in the matrix by a complex number and 
   * get a new matrix.
   */
  const matrix<DataType> operator*(const CD rhs) const;

  /**
   * Multiply each element in the matrix by a real number.
   */
  matrix<DataType> & operator*=(const double rhs);

  /**
   * Multiply each element in the matrix by a complex number.
   */
  matrix<DataType> & operator*=(const CD rhs);
 
  /**
   * Multiply two matrices and get a new matrix. Number of columns of the first
   * matrix must equals to the number of columns of the second matrix.
   */
  const matrix<DataType> operator*(const matrix<DataType>& rhs) const;

  /**
   * Multiply the matrix with a new matrix. Number of columns of the first
   * matrix must equals to the number of columns of the second matrix.
   */
  matrix<DataType>& operator*=(const matrix<DataType>& rhs);
 
  /**
   * Divide each element in the matrix by a constant to get a new matrix.
   */
  template<typename T> const matrix<DataType> operator/(const T rhs) const;

  /**
   * Divide each element in the matrix by a constant.
   */
  template<typename T> matrix<DataType>& operator/=(const T rhs);

  /*************************************
   *  matrix specific operations 
   ************************************/

  /**
   * Get the matrix's real part. 
   */
  matrix<double> real() const;

  /**
   * Get the matrix's imaginary part.  
   */
  matrix<double> imag() const;

  /**
   * Get the matrix's transposition. 
   */
  matrix<DataType> trans()const;

  /**
   * Get the matrix's hermitian.  
   */
  matrix<DataType> herm() const;

  /**
   * Get the matrix's conjugate.  
   */
  matrix<DataType> conj()const;

  /**
   * Get a sub-block of the matrix between r1th (containing) row and r2th 
   * (not containing) row, and between c1th (containing) column and c2th 
   * (not containing) column.  
   */
  matrix<DataType> block(size_t r1, size_t r2, size_t c1, size_t c2) const;

  /**
   * Replace a sub-block of the matrix with another matrix value. The
   * sub-block starts at (r-th row, c-th column) at its upper-left corner,
   * and is replace by values from mt's (0,0) position. 
   */
  matrix<DataType> &replace(size_t r, size_t c, const matrix<DataType> &mt);

  /**
   * Return the matrix's c-th column.
   */
  matrix<DataType> col(size_t c) const;

  /**
   * Return the matrix's r-th row.
   */
  matrix<DataType> row(size_t r) const;

  /**
   * Swap two cols.
   */
  matrix<DataType> &swapcol(size_t c1, size_t c2);

  /**
   * Swap two rows.
   */
  matrix<DataType> &swaprow(size_t r1, size_t r2);

  /**
   * Solve eigen system and put eigen values into e, corresponding 
   * left-eigen vectors into vl and corresponding right-eigen vectors
   * into vr. vl and vr will be square matrices.
   * 
   * If you are linking to basic blas library only, this function is not
   * implemented yet. 
   */
  void eigen(matrix<CD> &e, matrix<DataType> &vl, matrix<DataType> &vr);

  /**
   * Solve eigen system and put eigen values into e, corresponding 
   * left-eigen vectors into vl. vl will be square matrix.
   *
   * If you are linking to basic blas library only, this function is not
   * implemented yet. 
   */
  void reigen(matrix<CD> &e, matrix<DataType> &vr);

  /**
   * Solve eigen system and put eigen values into e, corresponding 
   * right-eigen vectors into vr. vr will be square matrix.
   *
   * If you are linking to basic blas library only, this function is not
   * implemented yet. 
   */
  void leigen(matrix<CD> &e, matrix<DataType> &vl);

/* private data */
protected:
  /**
   * Number of columns.
   */
  size_t m_ncol;

  /**
   * Number of rows.
   */
  size_t m_nrow;

};

typedef matrix<CD> cmatrix;
typedef matrix<double> dmatrix;

template<typename DataType> 
matrix<DataType>::matrix() : m_ncol(0), m_nrow(0) {
  m_temporary = false;
} 

template<typename DataType> 
matrix<DataType>::matrix(size_t nr, size_t nc) : m_ncol(nc), m_nrow(nr) {
  m_data = typename DataPtr<DataType>::Type(new DataArray<DataType>(m_ncol*m_nrow));
  m_temporary = false;
} 

template<typename DataType> 
template<typename T> 
matrix<DataType>::matrix(const T *data, size_t nr, size_t nc) 
  : m_ncol(nc), m_nrow(nr) {
  m_data =typename DataPtr<DataType>::Type(new DataArray<DataType>(data, nr*nc));
  m_temporary = false;
} 

template<typename DataType> 
template<typename T> 
matrix<DataType>::matrix(const matrix<T> &m) 
  : m_ncol(m.ncol()), m_nrow(m.nrow()) { 
  m_data =typename DataPtr<DataType>::Type(
    new DataArray<DataType>(m.getDataPtr(), m.nrow()*m.ncol())
  ); 
  m_temporary = false;
} 

template<typename DataType> 
matrix<DataType>::matrix(const matrix<DataType> &m) 
  : m_ncol(m.ncol()), m_nrow(m.nrow()) { 
  m_data =typename DataPtr<DataType>::Type(new DataArray<DataType>(*(m.m_data)));
  m_temporary = false;
}

template<typename DataType> 
matrix<DataType>& 
matrix<DataType>::operator= (const matrix<DataType> &rhs) {
  if (&rhs==this) return *this;
  if (m_ncol*m_nrow < rhs.ncol()*rhs.nrow()) 
    m_data =typename DataPtr<DataType>::Type(new DataArray<DataType>(*(rhs.m_data)));
  else 
    memcpy(m_data->m_data, rhs.getDataPtr(), m_ncol*m_nrow*sizeof(DataType));
  m_ncol = rhs.ncol();
  m_nrow = rhs.nrow();
  return *this;
} 

template<typename DataType> 
matrix<DataType>& 
matrix<DataType>::operator= (const DataType rhs) {
  size_t endi = m_ncol*m_nrow; 
  DataType *p = m_data->m_data;
  for (size_t i=0; i<endi; i++) *(p++)=rhs;
  return *this;
} 

template<typename DataType> 
bool matrix<DataType>::operator== (const matrix<DataType> &rhs) {
  if (&rhs==this) return true;
  if (m_nrow != rhs.nrow() || m_ncol != rhs.ncol()) return false;
  size_t endi = m_ncol*m_nrow; 
  DataType *p1 = m_data->m_data;
  DataType *p2 = rhs.getDataPtr();
  for (size_t i=0; i<endi; i++)
    if (*(p1++)!=*(p2++)) return false;
  return true;
} 

template<typename DataType> 
bool matrix<DataType>::operator!= (const matrix<DataType> &rhs) {
  return !(operator==(rhs));
}

/* interface to private data */
template<typename DataType> 
size_t matrix<DataType>::nrow() const {
  return m_nrow;
}

template<typename DataType> 
size_t matrix<DataType>::ncol() const {
  return m_ncol;
}

template<typename DataType> 
DataType * matrix<DataType>::getDataPtr() const {
  return m_data->m_data;
}

template<typename DataType> 
DataArray<DataType>& matrix<DataType>::getDataArray() const {
  return *m_data;
}

/* resizing/reshaping matrix */
template<typename DataType> 
void matrix<DataType>::resize(size_t nr, size_t nc) {
  if (m_nrow*m_ncol < nc * nr) {
    typename DataPtr<DataType>::Type 
      newData(new DataArray<DataType>(*m_data, nr*nc));
    m_data = newData;
  }
  m_ncol = nc;
  m_nrow = nr;
} 

template<typename DataType> 
inline void matrix<DataType>::clear() {
  if (m_data->m_data)
    memset(m_data->m_data, 0, m_data->m_size*sizeof(DataType));
}

/* operator overloading */
template<typename DataType> 
inline DataType& matrix<DataType>::operator()(size_t i, size_t j) {
  assert(i<m_nrow && j<m_ncol);
  return m_data->m_data[j*m_nrow + i];
}

template<typename DataType> 
inline DataType matrix<DataType>::operator()(size_t i, size_t j) const {
  assert(i<m_nrow && j<m_ncol);
  return m_data->m_data[j*m_nrow + i];
}

template<typename DataType> 
inline DataType& matrix<DataType>::operator()(size_t i) {
  return m_data->m_data[i];
}

template<typename DataType> 
inline DataType& matrix<DataType>::operator()(size_t i) const {
  return m_data->m_data[i];
}

template<typename DataType> 
inline DataType& matrix<DataType>::operator[](size_t i) {
  return m_data->m_data[i];
}

template<typename DataType> 
inline DataType& matrix<DataType>::operator[](size_t i) const {
  return m_data->m_data[i];
}

// arithmetic

template<typename DataType> 
template<typename T> 
const matrix<DataType> 
matrix<DataType>::operator+(const matrix<T>& rhs) const {
  assert(m_ncol == rhs.ncol() && m_nrow == rhs.nrow());
  matrix<DataType> m(m_nrow, m_ncol);
  size_t endi = m_nrow*m_ncol;
  DataType *p1, *p2, *p3;
  p1=m.getDataPtr();
  p2=getDataPtr();
  p3=rhs.getDataPtr();
  for (size_t i=0; i<endi; i++)
    *(p1++) = *(p2++) + *(p3++);
  return m;
}

template<typename DataType> 
template<typename T> 
matrix<DataType>& matrix<DataType>::operator+=(const matrix<T>& rhs) {
  assert(m_ncol == rhs.ncol() && m_nrow == rhs.nrow());
  size_t endi = m_nrow*m_ncol;
  DataType *p2, *p3;
  p2=getDataPtr();
  p3=rhs.getDataPtr();
  for (size_t i=0; i<endi; i++)
    *(p2++) += *(p3++);
  return *this;
}

template<typename DataType> 
template<typename T>
const matrix<DataType> 
matrix<DataType>::operator-(const matrix<T>& rhs) const {
  assert(m_ncol == rhs.ncol() && m_nrow == rhs.nrow());
  matrix<DataType> m(m_nrow, m_ncol);
  size_t endi = m_nrow*m_ncol;
  DataType *p1, *p2, *p3;
  p1=m.getDataPtr();
  p2=getDataPtr();
  p3=rhs.getDataPtr();
  for (size_t i=0; i<endi; i++)
    *(p1++) = *(p2++) - *(p3++);
  return m;
}

template<typename DataType> 
template<typename T> 
matrix<DataType>& matrix<DataType>::operator-=(const matrix<T>& rhs) {
  assert(m_ncol == rhs.ncol() && m_nrow == rhs.nrow());
  size_t endi = m_nrow*m_ncol;
  DataType *p2, *p3;
  p2=getDataPtr();
  p3=rhs.getDataPtr();
  for (size_t i=0; i<endi; i++)
    *(p2++) -= *(p3++);
  return *this;
}

template<typename DataType> 
const matrix<DataType> matrix<DataType>::operator*(const double rhs) const {
  matrix<DataType> m(m_nrow, m_ncol);
  size_t endi = m_nrow*m_ncol;
  DataType *p1, *p2;
  p1=m.getDataPtr();
  p2=getDataPtr();
  for (size_t i=0; i<endi; i++)
    *(p1++) = *(p2++)*rhs; 
  return m;
}

template<typename DataType> 
const matrix<DataType> matrix<DataType>::operator*(const CD rhs) const {
  matrix<DataType> m(m_nrow, m_ncol);
  size_t endi = m_nrow*m_ncol;
  DataType *p1, *p2;
  p1=m.getDataPtr();
  p2=getDataPtr();
  for (size_t i=0; i<endi; i++)
    *(p1++) = *(p2++)*rhs; 
  return m;
}

template<typename DataType> 
matrix<DataType> & matrix<DataType>::operator*=(const double rhs) {
  size_t endi = m_nrow*m_ncol;
  DataType *p2;
  p2=getDataPtr();
  for (size_t i=0; i<endi; i++)
    *(p2++) *= rhs; 
  return *this;
}

template<typename DataType> 
matrix<DataType> & matrix<DataType>::operator*=(const CD rhs) {
  size_t endi = m_nrow*m_ncol;
  DataType *p2;
  p2=getDataPtr();
  for (size_t i=0; i<endi; i++)
    *(p2++) *= rhs; 
  return *this;
}
 
template<typename DataType> 
const matrix<DataType> 
matrix<DataType>::operator*(const matrix<DataType>& rhs) const {
  assert (m_ncol == rhs.nrow());
  size_t m,n,k1,k2,lda,ldb;
  m = m_nrow;
  k1 = m_ncol;
  k2 = rhs.nrow();
  n = rhs.ncol();
  lda = m_nrow;
  ldb = rhs.nrow();

  matrix<DataType> ma;

  ma.resize(m,n);

  gemm<DataType>(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k1, 
    getDataPtr(), lda, rhs.getDataPtr(), ldb, ma.getDataPtr(), m);

  return ma;
} 

template<typename DataType> 
matrix<DataType>& matrix<DataType>::operator*=(const matrix<DataType>& rhs) {
  matrix<DataType> ma = (*this)*rhs; 
  (*this) = ma;
  return (*this);
}
 
template<typename DataType> 
template<typename T> 
const matrix<DataType> matrix<DataType>::operator/(const T rhs) const {
  matrix<DataType> m(m_nrow, m_ncol);
  size_t endi = m_nrow*m_ncol;
  DataType *p1, *p2;
  p1=m.getDataPtr();
  p2=getDataPtr();
  for (size_t i=0; i<endi; i++)
    *(p1++) = *(p2++)/rhs; 
  return m;
}

template<typename DataType> 
template<typename T> 
matrix<DataType>& matrix<DataType>::operator/=(const T rhs) {
  size_t endi = m_nrow*m_ncol;
  DataType *p2;
  p2=getDataPtr();
  for (size_t i=0; i<endi; i++)
    *(p2++) /= rhs; 
  return (*this);
}

/* matrix specific operations */
template<typename DataType> 
matrix<double> matrix<DataType>::real() const{
  matrix<double> m(m_nrow, m_ncol);
  size_t endi = m_nrow*m_ncol;
  double *p1;
  DataType *p2;
  p1=m.getDataPtr();
  p2=getDataPtr();
  for (size_t i=0; i<endi; i++) 
    *(p1++)=(p2++)->real;
  return m;
} 

template<typename DataType> 
matrix<double> matrix<DataType>::imag() const{
  matrix<double> m(m_nrow, m_ncol);
  size_t endi = m_nrow*m_ncol;
  double *p1;
  DataType *p2;
  p1=m.getDataPtr();
  p2=getDataPtr();
  for (size_t i=0; i<endi; i++) 
    *(p1++)=(p2++)->imag;
 return m;
}

template<typename DataType> 
matrix<DataType> matrix<DataType>::trans() const {
  matrix<DataType> m(m_ncol, m_nrow); 
  DataType *p0, *p1, *p2;
  p0=m.getDataPtr();
  p2=getDataPtr();
  for (size_t j=0; j<m_ncol; j++) { 
    p1 = p0 + j;
    for (size_t i=0; i<m_nrow; i++) { 
      *p1=*(p2++);
      p1+=m_ncol;
    }
  }
  return m;
} 

template<typename DataType> 
matrix<DataType> matrix<DataType>::herm() const {
  matrix<DataType> m(m_ncol, m_nrow); 
  DataType *p0, *p1, *p2;
  p0=m.getDataPtr();
  p2=getDataPtr();
  for (size_t j=0; j<m_ncol; j++) { 
    p1 = p0 + j;
    for (size_t i=0; i<m_nrow; i++) { 
      *p1=mlcpp::conj(*(p2++));
      p1+=m_ncol;
    }
  }
  return m;
} 

template<typename DataType> 
matrix<DataType> matrix<DataType>::conj() const {
  matrix<DataType> m(m_nrow, m_ncol); 
  size_t endi = m_nrow*m_ncol;
  DataType *p1, *p2;
  p1=m.getDataPtr();
  p2=getDataPtr();
  for (size_t i=0; i<endi; i++) 
      *(p1++)=mlcpp::conj(*(p2++));
  return m;
} 

template<typename DataType> 
matrix<DataType> matrix<DataType>::block(size_t r1, size_t r2, 
                                         size_t c1, size_t c2) const {
  assert(r1<=r2 && c1<=c2);
  matrix<DataType> m(r2-r1, c2-c1);
  DataType *p1, *p2;
  int n1, n2;
  n1 = m_nrow;
  n2 = n1-r2+r1;
  p2 = m.getDataPtr(); 
  p1 = m_data->m_data + r1+n1*c1;
  
  for (size_t j=c1; j<c2; j++) { 
    for (size_t i=r1; i<r2; i++) 
      *(p2++) = *(p1++);
    p1 += n2;
  }
  return m;
}

template<typename DataType> 
matrix<DataType> & matrix<DataType>::replace(size_t r, size_t c, 
                                            const matrix<DataType> &mt) {
  DataType *p1, *p2;
  int n1, n2, n3, n4;
  n1 = m_nrow;
  p2 = mt.getDataPtr(); 
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
matrix<DataType> matrix<DataType>::col(size_t c) const {
 return block(0,m_nrow, c, c+1);
}

template<typename DataType> 
matrix<DataType> matrix<DataType>::row(size_t r) const {
 return block(r, r+1, 0, m_ncol);
}

template<typename DataType> 
matrix<DataType> & matrix<DataType>::swapcol(size_t c1, size_t c2) {
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
matrix<DataType> & matrix<DataType>::swaprow(size_t r1, size_t r2) {
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
void matrix<DataType>::eigen(matrix<CD> &e, matrix<DataType> &vl, 
                             matrix<DataType> &vr) {
  assert (m_ncol == m_nrow);
  matrix<DataType> m(*this);

  vl.resize(m_ncol,m_ncol);
  vr.resize(m_ncol,m_ncol);
  e.resize(m_ncol, 1);
  
  char nl = &vl? 'V': 'N';
  char nr = &vr? 'V': 'N';

  geev<DataType>(nl, nr, m_ncol, m.getDataPtr(), m_ncol, 
                 e.getDataPtr(), vl.getDataPtr(), m_ncol, 
                 vr.getDataPtr(), m_ncol);
} 

template<typename DataType> 
void matrix<DataType>::reigen(matrix<CD> &e, matrix<DataType> &vr) {
  assert (m_ncol == m_nrow);
  matrix<DataType> m(*this);

  vr.resize(m_ncol,m_ncol);
  e.resize(m_ncol, 1);
  
  char nl = 'N';
  char nr = 'V';

  geev<DataType>(nl, nr, m_ncol, m.getDataPtr(), m_ncol, 
                 e.getDataPtr(), NULL, m_ncol, vr.getDataPtr(), m_ncol); 
} 

template<typename DataType> 
void matrix<DataType>::leigen(matrix<CD> &e, matrix<DataType> &vl) {
  assert (m_ncol == m_nrow);
  matrix<DataType> m(*this);

  vl.resize(m_ncol,m_ncol);
  e.resize(m_ncol, 1);
  
  char nl = 'V';
  char nr = 'N';

  geev<DataType>(nl, nr, m_ncol, m.getDataPtr(), m_ncol, 
                 e.getDataPtr(), vl.getDataPtr(), m_ncol, NULL, m_ncol); 
} 
/////////////////////////////////////////////////////////////////////////////

/** 
 * Return the real part of a matrix 
 */
matrix<double> real(const matrix<CD> &m) {
  return m.real();
} 

/** 
 * Return the imaginary part of a matrix 
 */
matrix<double> imag(const matrix<CD> &m) {
  return m.imag();
}

/** 
 * Return transpose of a matrix 
 */
template<typename DataType> 
matrix<DataType> trans(const matrix<DataType> &m) {
  return m.trans();
}

/** 
 * Return hermitian of a matrix 
 */
template<typename DataType> 
matrix<DataType> herm(const matrix<DataType> &m) {
  return m.herm();
}

/** 
 * Return conjugate of a matrix 
 */
template<typename DataType> 
matrix<DataType> conj(const matrix<DataType> &m) {
  return m.conj();
}

/** 
 * Multiply a double number and a matrix
 */
template<typename DataType> 
const matrix<DataType> operator*(const double lhs, 
                                 const matrix<DataType> &ma) {
  matrix<DataType> m(ma.nrow(), ma.ncol());
  size_t endi = ma.nrow()*ma.ncol();
  DataType *p1, *p2, *p3;
  p1=m.getDataPtr();
  p2=ma.getDataPtr();
  for (size_t i=0; i<endi; i++)
    *(p1++) = *(p2++)*lhs; 
  return m;
}

/** 
 * Multiply a complex number and a matrix
 */
template<typename DataType> 
const matrix<DataType> operator*(const CD lhs, 
                                 const matrix<DataType> &ma) {
  matrix<DataType> m(ma.nrow(), ma.ncol());
  size_t endi = ma.nrow()*ma.ncol();
  DataType *p1, *p2, *p3;
  p1=m.getDataPtr();
  p2=ma.getDataPtr();
  for (size_t i=0; i<endi; i++)
    *(p1++) = *(p2++)*lhs; 
  return m;
}

/* stream operator overload */
template<typename DataType> 
std::ostream& operator<< (std::ostream& out, 
                          const matrix<DataType> &m) {
  for (size_t i=0; i<m.nrow(); i++) {
    for (size_t j=0; j<m.ncol(); j++)
      out << m(i,j) << " ";
    if (i!=m.nrow()-1) out << std::endl;
  }
  return out; 
}

/**
 * Get print form of a matrix.
 */
template<typename DataType> 
std::string toString(const matrix<DataType> &m) {
  std::ostringstream out;
  out << m;
  return out.str(); 
}



/////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////
/** 
 * Identity matrix
 */
template<typename DataType> 
class idmatrix : public matrix<DataType> {
public:
  /**
   * Construct an identity matrix of n x n. 
   */
  idmatrix(size_t n);
};

typedef idmatrix<CD> cidmatrix;
typedef idmatrix<double> didmatrix;

template<typename DataType> 
idmatrix<DataType>::idmatrix(size_t n) {
  matrix<DataType>::m_data = 
    typename DataPtr<DataType>::Type(new DataArray<DataType>(n*n));
  matrix<DataType>::m_ncol = n;
  matrix<DataType>::m_nrow = n;
  
  size_t endi=n*n;
  DataType *p=matrix<DataType>::getDataPtr();
  for (size_t i=0; i<endi; i++) 
    *(p++) = 0;

  p=matrix<DataType>::getDataPtr();
  for (size_t i=0; i<n; i++) 
    *(p+i*(n+1)) = 1.;

}
/////////////////////////////////////////////////////////////////////////////


}
#endif
