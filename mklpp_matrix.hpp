#ifndef MKLPP_MATRIX_H
#define MKLPP_MATRIX_H
#include <mklpp_complex.hpp>
#include <mklpp_wrapper.hpp>
#include <mklpp_dataarray.hpp>

namespace mklpp {

/////////////////////////////////////////////////////////////////////////////
/* matrix class */
template<typename DataType> 
class matrix {
public:
  size_t nCol;
  size_t nRow;

  ~matrix() {}

  matrix();

  matrix(size_t nr, size_t nc);
 
  template<typename T> 
  matrix(const T *data, size_t nr, size_t nc); 

  template<typename T> 
  matrix(const matrix<T> &m);
 
  matrix(const matrix<DataType> &m); 

  const matrix<DataType>& operator= (const matrix<DataType> &rhs);
 
  const matrix<DataType>& operator= (const DataType rhs);
 
  bool operator== (const matrix<DataType> &rhs); 

/* interface to private data */
  DataType * getDataPtr() const; 
  DataArray<DataType>& getDataArray() const; 

  size_t getDataSize() const; 

  CBLAS_TRANSPOSE getTranspose() const;
 
/* resizing/reshaping matrix */
  void resize(size_t nr, size_t nc); 

/* indexing */
  inline DataType& operator()(size_t i, size_t j); 

  inline DataType operator()(size_t i, size_t j) const; 

  inline DataType& operator()(size_t i);
 
  inline DataType& operator()(size_t i) const;

  inline DataType& operator[](size_t i);
 
  inline DataType& operator[](size_t i) const; 

/* arithmetic */
  template<typename T> 
  const matrix<DataType> operator+(const matrix<T>& rhs) const;
 
  template<typename T> matrix<DataType>& operator+=(const matrix<T>& rhs);
 
  template<typename T> const matrix<DataType> operator-(const matrix<T>& rhs) const;

  template<typename T> matrix<DataType>& operator-=(const matrix<T>& rhs);

  const matrix<DataType> operator*(const double rhs);

  const matrix<DataType> operator*(const CD rhs);

  matrix<DataType> & operator*=(const double lhs);

  matrix<DataType> & operator*=(const CD lhs);
 
  const matrix<DataType> operator*(const matrix<DataType>& rhs) const;

  matrix<DataType>& operator*=(const matrix<DataType>& rhs);
 
  template<typename T> const matrix<DataType> operator/(const T rhs);

  template<typename T> matrix<DataType>& operator/=(const T rhs);

/* matrix specific operations */
  matrix<DataType> &trans();

  matrix<DataType> &herm();

  matrix<DataType> &conj();

  matrix<DataType> block(size_t r1, size_t r2, size_t c1, size_t c2);

  matrix<DataType> &insert(size_t r, size_t c, matrix<DataType> mt);

  matrix<DataType> col(size_t c);

  matrix<DataType> row(size_t r);

/* solvers */
  void eigen(matrix<CD> &e, matrix<DataType> &vl, matrix<DataType> &vr);

  void reigen(matrix<CD> &e, matrix<DataType> &vr);

  void leigen(matrix<CD> &e, matrix<DataType> &vl);

/* private data */
protected:
  typename DataPtr<DataType>::Type _Data; 
  size_t _DataSize;
  CBLAS_TRANSPOSE _Transpose;
};

typedef matrix<CD> cmatrix;
typedef matrix<double> dmatrix;

template<typename DataType> 
matrix<DataType>::matrix() : nCol(0), nRow(0) {
  _DataSize = 0;
  _Transpose = CblasNoTrans;
} 

template<typename DataType> 
matrix<DataType>::matrix(size_t nr, size_t nc) : nCol(nc), nRow(nr) {
    _Data = DataPtr<DataType>::Type(new DataArray<DataType>(nCol*nRow));
    _DataSize = nr*nc;
    _Transpose = CblasNoTrans;
  } 

template<typename DataType> 
template<typename T> 
matrix<DataType>::matrix(const T *data, size_t nr, size_t nc) 
  : nCol(nc), nRow(nr) {
  _Data = DataPtr<DataType>::Type(new DataArray<DataType>(data, nr*nc));
  _DataSize = nr*nc;
  _Transpose = CblasNoTrans;
} 

template<typename DataType> 
template<typename T> 
matrix<DataType>::matrix(const matrix<T> &m) 
  : nCol(m.nCol), nRow(m.nRow), 
  _DataSize(m.nRow*m.nCol), _Transpose(m.getTranspose()){
  _Data = DataPtr<DataType>::Type(
    new DataArray<DataType>(m.getDataPtr(), m.nRow*m.nCol)
  ); 
} 

template<typename DataType> 
matrix<DataType>::matrix(const matrix<DataType> &m) 
  : nCol(m.nCol), nRow(m.nRow), 
  _DataSize(m._DataSize), _Transpose(m._Transpose){
  _Data = DataPtr<DataType>::Type(new DataArray<DataType>(*(m._Data)));
}

template<typename DataType> 
const matrix<DataType>& 
matrix<DataType>::operator= (const matrix<DataType> &rhs) {
  if (&rhs==this) return *this;
  nCol = rhs.nCol;
  nRow = rhs.nRow;
  if (nCol*nRow > _DataSize) {
    _Data = rhs._Data;
    _DataSize = nCol*nRow;
  }
  _Transpose = rhs._Transpose;
  return *this;
} 

template<typename DataType> 
const matrix<DataType>& 
matrix<DataType>::operator= (const DataType rhs) {
  for (size_t i=0; i<_DataSize; i++) (*_Data)[i]=rhs;
  return *this;
} 

template<typename DataType> 
bool matrix<DataType>::operator== (const matrix<DataType> &rhs) {
  if (nRow != rhs.nRow || nCol != rhs.nCol) return false;
  for (size_t i=0; i<nRow; i++)
    for (size_t j=0; j<nCol; j++)
      if ((*this)(i,j)!=rhs(i,j)) return false;
  return true;
} 

/* interface to private data */
template<typename DataType> 
DataType * matrix<DataType>::getDataPtr() const {
  return (*_Data)._Data;
}

template<typename DataType> 
DataArray<DataType>& matrix<DataType>::getDataArray() const {
  return *_Data;
}

template<typename DataType> 
size_t matrix<DataType>::getDataSize() const {
  return _DataSize;
}

template<typename DataType> 
CBLAS_TRANSPOSE matrix<DataType>::getTranspose() const {
  return _Transpose;
}

/* resizing/reshaping matrix */
template<typename DataType> 
void matrix<DataType>::resize(size_t nr, size_t nc) {
  if (_DataSize < nc * nr) {
    typename DataPtr<DataType>::Type 
      newData(new DataArray<DataType>(*_Data, nr*nc));
    _Data = newData;
    _DataSize = nc*nr;
  }
  nCol = nc;
  nRow = nr;
} 

/* operator overloading */
template<typename DataType> 
inline DataType& matrix<DataType>::operator()(size_t i, size_t j) {
  assert(i<nRow && j<nCol);
  if (_Transpose == CblasNoTrans) return (*_Data)[j*nRow + i];
  else return (*_Data)[i*nCol + j];
}

template<typename DataType> 
inline DataType matrix<DataType>::operator()(size_t i, size_t j) const {
  assert(i<nRow && j<nCol);
  if (_Transpose == CblasNoTrans) return (*_Data)[j*nRow + i];
  else return (*_Data)[i*nCol + j];
}

template<typename DataType> 
inline DataType& matrix<DataType>::operator()(size_t i) {
  assert((nCol == 1 && i<nRow) || (nRow == 1 && i<nCol));
  return (*_Data)[i];
}

template<typename DataType> 
inline DataType& matrix<DataType>::operator()(size_t i) const {
  assert((nCol == 1 && i<nRow) || (nRow == 1 && i<nCol));
  return (*_Data)[i];
}

template<typename DataType> 
inline DataType& matrix<DataType>::operator[](size_t i) {
  assert((nCol == 1 && i<nRow) || (nRow == 1 && i<nCol));
  return (*_Data)[i];
}

template<typename DataType> 
inline DataType& matrix<DataType>::operator[](size_t i) const {
  assert((nCol == 1 && i<nRow) || (nRow == 1 && i<nCol));
  return (*_Data)[i];
}

// arithmetic

template<typename DataType> 
template<typename T> 
const matrix<DataType> 
matrix<DataType>::operator+(const matrix<T>& rhs) const {
  assert(nCol == rhs.nCol && nRow == rhs.nRow);
  matrix<DataType> m(nRow, nCol);
  for (size_t i=0; i<nRow; i++)
    for (size_t j=0; j<nCol; j++)
      m(i,j) = (*this)(i,j)+rhs(i,j); 
  return m;
}

template<typename DataType> 
template<typename T> 
matrix<DataType>& matrix<DataType>::operator+=(const matrix<T>& rhs) {
  assert(nCol == rhs.nCol && nRow == rhs.nRow);
  for (size_t i=0; i<nRow; i++)
    for (size_t j=0; j<nCol; j++)
      (*this)(i,j) += rhs(i,j); 
  return *this;
}

template<typename DataType> 
template<typename T>
const matrix<DataType> 
matrix<DataType>::operator-(const matrix<T>& rhs) const {
  assert(nCol == rhs.nCol && nRow == rhs.nRow);
  matrix<DataType> m(nRow, nCol);
  for (size_t i=0; i<nRow; i++)
    for (size_t j=0; j<nCol; j++)
      m(i,j) = (*this)(i,j)-rhs(i,j); 
  return m;
}

template<typename DataType> 
template<typename T> 
matrix<DataType>& matrix<DataType>::operator-=(const matrix<T>& rhs) {
  assert(nCol == rhs.nCol && nRow == rhs.nRow);
  for (size_t i=0; i<nRow; i++)
    for (size_t j=0; j<nCol; j++)
      (*this)(i,j) -= rhs(i,j); 
  return *this;
}

template<typename DataType> 
const matrix<DataType> matrix<DataType>::operator*(const double rhs) {
  matrix<DataType> m(nRow, nCol);
  for (size_t i=0; i<nRow; i++)
    for (size_t j=0; j<nCol; j++)
      m(i,j) = (*this)(i,j)*rhs; 
  return m;
}

template<typename DataType> 
const matrix<DataType> matrix<DataType>::operator*(const CD rhs) {
  matrix<DataType> m(nRow, nCol);
  for (size_t i=0; i<nRow; i++)
    for (size_t j=0; j<nCol; j++)
      m(i,j) = (*this)(i,j)*rhs; 
  return m;
}

template<typename DataType> 
matrix<DataType> & matrix<DataType>::operator*=(const double lhs) {
  for (size_t i=0; i<ma.nRow; i++)
    for (size_t j=0; j<ma.nCol; j++)
    (*this)(i,j) *= lhs; 
  return *this;
}

template<typename DataType> 
matrix<DataType> & matrix<DataType>::operator*=(const CD lhs) {
  for (size_t i=0; i<ma.nRow; i++)
    for (size_t j=0; j<ma.nCol; j++)
    (*this)(i,j) *= lhs; 
  return *this;
}
 
template<typename DataType> 
const matrix<DataType> 
matrix<DataType>::operator*(const matrix<DataType>& rhs) const {
  assert (nCol == rhs.nRow);
  size_t m,n,k1,k2,lda,ldb;
  m = nRow;
  k1 = nCol;
  k2 = rhs.nRow;
  n = rhs.nCol;
  lda = (_Transpose==CblasTrans)? nCol: nRow;
  ldb = (rhs._Transpose==CblasTrans)? rhs.nCol: rhs.nRow;

  matrix<DataType> ma;

  ma.resize(m,n);

  gemm<DataType>(CblasColMajor, _Transpose, rhs._Transpose, m, n, k1, 
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
const matrix<DataType> matrix<DataType>::operator/(const T rhs) {
  matrix<DataType> m(nRow, nCol);
  for (size_t i=0; i<nRow; i++)
    for (size_t j=0; j<nCol; j++)
      m(i,j) = (*this)(i,j)/rhs; 
  return m;
}

template<typename DataType> 
template<typename T> 
matrix<DataType>& matrix<DataType>::operator/=(const T rhs) {
  for (size_t i=0; i<nRow; i++)
    for (size_t j=0; j<nCol; j++)
      (*this)(i,j) /= rhs; 
  return (*this);
}

/* matrix specific operations */
template<typename DataType> 
matrix<DataType> & matrix<DataType>::trans() {
  _Transpose = (_Transpose == CblasTrans)? CblasNoTrans: CblasTrans;
 resize(nCol, nRow);
 return *this;
} 

template<typename DataType> 
matrix<DataType> & matrix<DataType>::herm() {
 trans();
 for (size_t i=0; i<nRow; i++) 
   for (size_t j=0; j<nCol; j++) 
     (*this)(i,j)=mklpp::conj((*this)(i,j));
 return *this;
} 

template<typename DataType> 
matrix<DataType> & matrix<DataType>::conj() {
 for (size_t i=0; i<nRow; i++) 
   for (size_t j=0; j<nCol; j++) 
     (*this)(i,j)=mklpp::conj((*this)(i,j));
 return *this;
} 

template<typename DataType> 
matrix<DataType> matrix<DataType>::block(size_t r1, size_t r2, 
                                         size_t c1, size_t c2) {
  assert(r1<=r2 && c1<=c2);
  matrix<DataType> m(r2-r1, c2-c1);
  for (size_t i=0; i<m.nRow; i++) 
    for (size_t j=0; j<m.nCol; j++) 
      m(i,j) = (*this)(i+r1, j+c1);
  return m;
}

template<typename DataType> 
matrix<DataType> & matrix<DataType>::insert(size_t r, size_t c, 
                                            matrix<DataType> mt) {
  for (size_t i=0; i<mt.nRow; i++) 
    for (size_t j=0; j<mt.nCol; j++) 
      (*this)(i+r, j+c) = mt(i,j);
  return *this;
}

template<typename DataType> 
matrix<DataType> matrix<DataType>::col(size_t c) {
 return block(0,nRow, c, c+1);
}

template<typename DataType> 
matrix<DataType> matrix<DataType>::row(size_t r) {
 return block(r, r+1, 0, nCol);
}

/* solvers */
template<typename DataType> 
void matrix<DataType>::eigen(matrix<CD> &e, matrix<DataType> &vl, 
                             matrix<DataType> &vr) {
  assert (nCol == nRow && 
          (_Transpose == CblasNoTrans || _Transpose == CblasTrans));
  matrix<DataType> m(*this);
  if ( _Transpose == CblasTrans) {
    m.resize(nCol, nCol);
    for (size_t i=0; i<nRow; i++) 
      for (size_t j=0; j<nCol; j++) 
        m(i,j)= (*this)(j,i);
  }    

  vl.resize(nCol,nCol);
  vr.resize(nCol,nCol);
  e.resize(nCol, 1);
  
  char nl = &vl? 'V': 'N';
  char nr = &vr? 'V': 'N';

  geev<DataType>(LAPACK_COL_MAJOR, nl, nr, nCol, m.getDataPtr(), nCol, 
                 e.getDataPtr(), vl.getDataPtr(), nCol, 
                 vr.getDataPtr(), nCol);
} 

template<typename DataType> 
void matrix<DataType>::reigen(matrix<CD> &e, matrix<DataType> &vr) {
  assert (nCol == nRow && (_Transpose == CblasNoTrans || 
                           _Transpose == CblasTrans));
  matrix<DataType> m(*this);
  if ( _Transpose == CblasTrans) {
    m.resize(nCol, nCol);
    for (size_t i=0; i<nRow; i++) 
      for (size_t j=0; j<nCol; j++) 
        m(i,j) = (*this)(j,i);
  }    

  vr.resize(nCol,nCol);
  e.resize(nCol, 1);
  
  char nl = 'N';
  char nr = 'V';

  geev<DataType>(LAPACK_COL_MAJOR, nl, nr, nCol, m.getDataPtr(), nCol, 
                 e.getDataPtr(), NULL, nCol, vr.getDataPtr(), nCol); 
} 

template<typename DataType> 
void matrix<DataType>::leigen(matrix<CD> &e, matrix<DataType> &vl) {
  assert (nCol == nRow && (_Transpose == CblasNoTrans || 
                           _Transpose == CblasTrans));
  matrix<DataType> m(*this);
  if ( _Transpose == CblasTrans) {
    m.resize(nCol, nCol);
    for (size_t i=0; i<nRow; i++) 
      for (size_t j=0; j<nCol; j++) 
        m(i,j) = (*this)(j,i);
  }    

  vl.resize(nCol,nCol);
  e.resize(nCol, 1);
  
  char nl = 'V';
  char nr = 'N';

  geev<DataType>(LAPACK_COL_MAJOR, nl, nr, nCol, m.getDataPtr(), nCol, 
                 e.getDataPtr(), vl.getDataPtr(), nCol, NULL, nCol); 
} 
/////////////////////////////////////////////////////////////////////////////

matrix<double> real(const matrix<CD> &ma) {
 matrix<double> m(ma.nRow, ma.nCol);
 for (size_t i=0; i<ma.nRow; i++) 
   for (size_t j=0; j<ma.nCol; j++) 
     m(i,j)=ma(i,j).real;
 return m;
} 

matrix<double> imag(const matrix<CD> &ma) {
 matrix<double> m(ma.nRow, ma.nCol);
 for (size_t i=0; i<ma.nRow; i++) 
   for (size_t j=0; j<ma.nCol; j++) 
     m(i,j)=ma(i,j).imag;
 return m;
}

template<typename DataType> 
matrix<DataType> trans(const matrix<DataType> &m) {
  matrix<DataType> m1 = m;
  m.trans();
  return m1;
}

template<typename DataType> 
matrix<DataType> herm(const matrix<DataType> &m) {
  matrix<DataType> m1 = m;
  m1.herm();
  return m1;
}

template<typename DataType> 
matrix<DataType> conj(const matrix<DataType> &m) {
  matrix<DataType> m1 = m;
  m.conj();
  return m1;
}

/* arithmetic */
template<typename DataType> 
const matrix<DataType> operator*(const double lhs, 
                                 const matrix<DataType> &ma) {
  matrix<DataType> m(ma.nRow, ma.nCol);
  for (size_t i=0; i<ma.nRow; i++)
    for (size_t j=0; j<ma.nCol; j++)
      m(i,j) = ma(i,j)*lhs; 
  return m;
}

template<typename DataType> 
const matrix<DataType> operator*(const CD lhs, 
                                 const matrix<DataType> &ma) {
  matrix<DataType> m(ma.nRow, ma.nCol);
  for (size_t i=0; i<ma.nRow; i++)
    for (size_t j=0; j<ma.nCol; j++)
      m(i,j) = ma(i,j)*lhs; 
  return m;
}

/* stream operator overload */
template<typename DataType> 
std::ostream& operator<< (std::ostream& out, 
                          const matrix<DataType> &m) {
  for (size_t i=0; i<m.nRow; i++) {
    for (size_t j=0; j<m.nCol; j++)
      out << m(i,j) << " ";
    if (i!=m.nRow-1) out << std::endl;
  }
  return out; 
}

/////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////
/* idmatrix class
   identity matrix */
template<typename DataType> 
class idmatrix : public matrix<DataType> {
public:
  idmatrix(size_t n);
};

typedef idmatrix<CD> cidmatrix;
typedef idmatrix<double> didmatrix;

template<typename DataType> 
idmatrix<DataType>::idmatrix(size_t n) {
  matrix<DataType>::_Data = 
    DataPtr<DataType>::Type(new DataArray<DataType>(n*n));
  matrix<DataType>::nCol = n;
  matrix<DataType>::nRow = n;
  matrix<DataType>::_DataSize = n*n;
  matrix<DataType>::_Transpose = CblasNoTrans;
  
  for (size_t i=0; i<n; i++) 
    for (size_t j=0; j<n; j++)
      if (i==j) (*this)(i,j)=1.;
      else (*this)(i,j)=0.;
}
/////////////////////////////////////////////////////////////////////////////


}
#endif
