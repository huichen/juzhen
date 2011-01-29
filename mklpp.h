#ifndef MKLPP_H
#include <iostream>
#include <mkl_cblas.h>
#include <mkl.h>
#include <string.h>
#include <assert.h>
#include <complex>
#include <cmath>
#include <memory>

namespace mklpp {
using namespace std;

#define CD complex<double>
#define cmatrix matrix<complex<double> >
#define cidmatrix idmatrix<complex<double> >

/* complex constant */
complex<double> C1 = complex<double>(1,0);
complex<double> Ci = complex<double>(0,1);

/* MKL function wrappers implemented with templates*/
template<typename T> void gemm(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N, const MKL_INT K, const T alpha, const T *A, const MKL_INT lda, const T *B, const MKL_INT ldb, const T beta, T *c, const MKL_INT ldc) { 
  assert(0); // always fails
};

template<> void gemm<double>(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N, const MKL_INT K, const double alpha, const double *A, const MKL_INT lda, const double *B, const MKL_INT ldb, const double beta, double *c, const MKL_INT ldc) { 
  cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, c, ldc); 
};

template<> void gemm<CD >(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N, const MKL_INT K, const CD  alpha, const CD  *A, const MKL_INT lda, const CD  *B, const MKL_INT ldb, const CD  beta, CD  *c, const MKL_INT ldc) { 
  cblas_zgemm(Order, TransA, TransB, M, N, K, &alpha, A, lda, B, ldb, &beta, c, ldc); 
};

template<typename DataType> class DataArray {
public:
  DataArray(size_t s) {
    _Data = (DataType *) malloc(s*sizeof(DataType));
    assert (_Data);
    size = s;
  }

  DataArray(const DataArray<DataType> &da) {
    if (&da) {
      _Data = (DataType *) malloc(da.size*sizeof(DataType));
      assert (_Data);
      if (da._Data) memcpy(_Data, da._Data, da.size*sizeof(DataType));
      size = da.size;
    } else {
      size = 0;
      _Data = NULL;
    }
  }

  DataArray(const DataArray<DataType> &da, const size_t s) {
    _Data = (DataType *) malloc(s*sizeof(DataType));
    assert (_Data);
    if (&da) {
      size = (s < da.size? s : da.size);
      memcpy(_Data, da._Data, size*sizeof(DataType));
    }
  }

  DataArray(const DataType *da, const size_t s) {
    _Data = (DataType *) malloc(s*sizeof(DataType));
    assert (_Data);
    if (da) memcpy(_Data, da, s*sizeof(DataType));
    size = s;
  }
   
  ~DataArray() {
    if (_Data) free(_Data);
  } 

  DataType& operator[](const size_t i) {
    return _Data[i];
  }
 
  DataType * _Data;
  size_t size;
};

template<typename DataType>
class DataPtr{
public:
  typedef auto_ptr<DataArray<DataType> > Type;
};

/* matrix class */
template<typename DataType> class matrix {
public:
  size_t nCol;
  size_t nRow;

/* deconstruction, construction and assignment */
  ~matrix() {
  }

  matrix() : nCol(0), nRow(0) {
    _DataSize = 0;
    _Transpose = CblasNoTrans;
  } 

  matrix(size_t nr, size_t nc) : nCol(nc), nRow(nr) {
    _Data = DataPtr<DataType>::Type(new DataArray<DataType>(nCol*nRow));
    _DataSize = nr*nc;
    _Transpose = CblasNoTrans;
  } 

  matrix(const DataType *data, size_t nr, size_t nc) : nCol(nc), nRow(nr) {
    _Data = DataPtr<DataType>::Type(new DataArray<DataType>(data, nr*nc));
    _DataSize = nr*nc;
    _Transpose = CblasNoTrans;
  } 

  matrix(const matrix<DataType> &m) : nCol(m.nCol), nRow(m.nRow), _DataSize(m._DataSize), _Transpose(m._Transpose){
    _Data = DataPtr<DataType>::Type(new DataArray<DataType>(*(m._Data))); 
  } 

  const matrix<DataType>& operator= (const matrix<DataType> &rhs) {
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

  const matrix<DataType>& operator= (const DataType rhs) {
    for (size_t i=0; i<_DataSize; i++) (*_Data)[i]=rhs;
    return *this;
  } 

/* interface to private data */
  DataType * getDataPtr() const {
    return (*_Data)._Data;
  }

  size_t getDataSize() {
    return _DataSize;
  }

/* resizing/reshaping matrix */
  void resize(size_t nr, size_t nc) {
    if (_DataSize < nc * nr) {
      typename DataPtr<DataType>::Type newData(new DataArray<DataType>(*_Data, nr*nc));
      _Data = newData;
      _DataSize = nc*nr;
    }
    nCol = nc;
    nRow = nr;
  } 

/* operator overloading */
  inline DataType& operator()(size_t i, size_t j) {
    assert(i<nRow && j<nCol);
    if (_Transpose == CblasNoTrans) return (*_Data)[j*nRow + i];
    else return (*_Data)[i*nCol + j];
  }

  inline DataType operator()(size_t i, size_t j) const {
    assert(i<nRow && j<nCol);
    if (_Transpose == CblasNoTrans) return (*_Data)[j*nRow + i];
    else return (*_Data)[i*nCol + j];
  }

  // arithmetic
  const matrix<DataType> operator*(const matrix<DataType>& rhs) const {
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
    DataType alpha, beta;
    alpha = 1;
    beta = 0;
    
    gemm<DataType>(CblasColMajor, _Transpose, rhs._Transpose, m, n, k1, 
      alpha, (*_Data)._Data, lda, rhs.getDataPtr(), ldb, beta, ma.getDataPtr(), m);

    return ma;
  } 

  template<typename T> const matrix<DataType> operator+(const matrix<T>& rhs) {
    assert(nCol == rhs.nCol && nRow == rhs.nRow);
    matrix<DataType> m(nRow, nCol);
    for (size_t i=0; i<nRow; i++)
      for (size_t j=0; j<nCol; j++)
        m(i,j) = (*this)(i,j)+rhs(i,j); 
    return m;
  }

  template<typename T> matrix<DataType>& operator+=(const matrix<T>& rhs) {
    assert(nCol == rhs.nCol && nRow == rhs.nRow);
    for (size_t i=0; i<nRow; i++)
      for (size_t j=0; j<nCol; j++)
        (*this)(i,j) += rhs(i,j); 
    return *this;
  }

  template<typename T> const matrix<DataType> operator-(const matrix<T>& rhs) {
    assert(nCol == rhs.nCol && nRow == rhs.nRow);
    matrix<DataType> m(nRow, nCol);
    for (size_t i=0; i<nRow; i++)
      for (size_t j=0; j<nCol; j++)
        m(i,j) = (*this)(i,j)-rhs(i,j); 
    return m;
  }

  template<typename T> matrix<DataType>& operator-=(const matrix<T>& rhs) {
    assert(nCol == rhs.nCol && nRow == rhs.nRow);
    for (size_t i=0; i<nRow; i++)
      for (size_t j=0; j<nCol; j++)
        (*this)(i,j) -= rhs(i,j); 
    return *this;
  }

  template<typename T> const matrix<DataType> operator/(const T rhs) {
    matrix<DataType> m(nRow, nCol);
    for (size_t i=0; i<nRow; i++)
      for (size_t j=0; j<nCol; j++)
        m(i,j) = (*this)(i,j)/rhs; 
    return m;
  }

  template<typename T> matrix<DataType>& operator/=(const T rhs) {
    for (size_t i=0; i<nRow; i++)
      for (size_t j=0; j<nCol; j++)
        (*this)(i,j) /= rhs; 
    return (*this);
  }

/* matrix specific operations */
  matrix<DataType> &trans() {
    _Transpose = (_Transpose == CblasTrans)? CblasNoTrans: CblasTrans;
   resize(nCol, nRow);
   return *this;
  } 

  matrix<DataType> &herm() {
   trans();
   for (size_t i=0; i<nRow; i++) 
     for (size_t j=0; j<nCol; j++) 
       (*this)(i,j)=std::conj((*this)(i,j));
   return *this;
  } 

  matrix<DataType> &conj() {
   for (size_t i=0; i<nRow; i++) 
     for (size_t j=0; j<nCol; j++) 
       (*this)(i,j)=std::conj((*this)(i,j));
   return *this;
  } 

  matrix<DataType> block(size_t r1, size_t r2, size_t c1, size_t c2) {
    assert(r1<=r2 && c1<=c2);
    matrix<DataType> m(r2-r1, c2-c1);
    for (size_t i=0; i<m.nRow; i++) 
      for (size_t j=0; j<m.nCol; j++) 
        m(i,j) = (*this)(i+r1, j+c1);
    return m;
  }

  matrix<DataType> &insert(size_t r, size_t c, matrix<DataType> mt) {
    for (size_t i=0; i<mt.nRow; i++) 
      for (size_t j=0; j<mt.nCol; j++) 
        (*this)(i+r, j+c) = mt(i,j);
    return *this;
  }

  matrix<DataType> col(size_t c) {
   return block(0,nRow, c, c+1);
  }

  matrix<DataType> row(size_t r) {
   return block(r, r+1, 0, nCol);
  }

/* private data */
protected:
  typename DataPtr<DataType>::Type _Data; 
  size_t _DataSize;
  CBLAS_TRANSPOSE _Transpose;
};

/* matrix derived classes */
template<typename DataType> class idmatrix : public matrix<DataType> {
// identity matrix
public:
  idmatrix(size_t n) {
    matrix<DataType>::_Data = DataPtr<DataType>::Type(new DataArray<DataType>(n*n));
    matrix<DataType>::nCol = n;
    matrix<DataType>::nRow = n;
    matrix<DataType>::_DataSize = n*n;
    
    for (size_t i=0; i<n; i++) 
      for (size_t j=0; j<n; j++)
        if (i==j) (*this)(i,j)=1.;
        else (*this)(i,j)=0;
  }
};

/* matrix operation functions */

template<typename DataType> matrix<DataType> trans(const matrix<DataType> &m) {
  matrix<DataType> m1 = m;
  m.trans();
  return m1;
}

template<typename DataType> matrix<DataType> herm(const matrix<DataType> &m) {
  matrix<DataType> m1 = m;
  m.herm();
  return m1;
}

template<typename DataType> matrix<DataType> conj(const matrix<DataType> &m) {
  matrix<DataType> m1 = m;
  m.conj();
  return m1;
}

/* arithmetic */
template<typename DataType, typename T> const matrix<DataType> operator*(const T lhs, const matrix<DataType> &ma) {
  matrix<DataType> m(ma.nRow, ma.nCol);
  for (size_t i=0; i<ma.nRow; i++)
    for (size_t j=0; j<ma.nCol; j++)
      m(i,j) = ma(i,j)*lhs; 
  return m;
}

template<typename DataType> const matrix<DataType> operator*(const matrix<DataType>& lhs, const matrix<DataType>& ma) {
  return ma.operator*(lhs);
}

/* stream operator overload */
template<typename DataType> ostream& operator<< (ostream& out, const matrix<DataType> &m) {
  for (size_t i=0; i<m.nRow; i++) {
    for (size_t j=0; j<m.nCol; j++)
      out << m(i,j) << " ";
    if (i!=m.nRow-1) out << endl;
  }
  return out; 
}
 
}
#endif
