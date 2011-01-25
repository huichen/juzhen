#ifndef MKLPP_H
#include <iostream>
#include <mkl_types.h>
#include <mkl_cblas.h>
#include <mkl.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <complex>
#include <cmath>

namespace MAT {
using namespace std;

#define CD complex<double>
#define CMatrix Matrix<complex<double> >
#define CDMatrix DMatrix<complex<double> >

/* complex constant */
complex<double> C1 = complex<double>(1,0);
complex<double> Ci = complex<double>(0,1);

/* MKL function wrappers implemented with templates*/
template<typename T> void gemm(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N, const MKL_INT K, const T alpha, const T *A, const MKL_INT lda, const T *B, const MKL_INT ldb, const T beta, T *c, const MKL_INT ldc) { };

template<> void gemm<double>(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N, const MKL_INT K, const double alpha, const double *A, const MKL_INT lda, const double *B, const MKL_INT ldb, const double beta, double *c, const MKL_INT ldc) { 
  cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, c, ldc); 
};

template<> void gemm<complex<double> >(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N, const MKL_INT K, const complex<double>  alpha, const complex<double>  *A, const MKL_INT lda, const complex<double>  *B, const MKL_INT ldb, const complex<double>  beta, complex<double>  *c, const MKL_INT ldc) { 
  cblas_zgemm(Order, TransA, TransB, M, N, K, &alpha, A, lda, B, ldb, &beta, c, ldc); 
};

/* Matrix class */
template<typename DataType> class Matrix {
public:
  size_t nCol;
  size_t nRow;

/* deconstruction, construction and assignment */
  ~Matrix() {
    if(_Data) free(_Data);
  }

  Matrix() : nCol(0), nRow(0) {
    _Data = NULL;  
    _DataSize = 0;
    _Transpose = CblasNoTrans;
  } 

  Matrix(size_t nr, size_t nc) : nCol(nc), nRow(nr) {
    _Data = (DataType *) malloc(nCol*nRow*sizeof(DataType));
    assert(_Data);
    _DataSize = nr*nc;
    _Transpose = CblasNoTrans;
  } 

  Matrix(const DataType *data, size_t nr, size_t nc) : nCol(nc), nRow(nr) {
    nCol = nc;
    nRow = nr;
    _Data = (DataType *) malloc(nCol*nRow*sizeof(DataType));
    assert(_Data);
    memcpy(_Data, data, nCol*nRow*sizeof(DataType));    
    _DataSize = nCol * nRow;
    _Transpose = CblasNoTrans;
  } 

  Matrix(const Matrix<DataType> &m) {
    nCol = m.nCol;
    nRow = m.nRow;
    _Data = (DataType *) malloc(nCol*nRow*sizeof(DataType));
    memcpy(_Data, m.getDataPtr(), nCol*nRow*sizeof(DataType));    
    _DataSize = nCol*nRow;
    _Transpose = m._Transpose;
  } 

  const Matrix<DataType>& operator= (const Matrix<DataType> &rhs) {
    if (&rhs==this) return *this;
    nCol = rhs.nCol;
    nRow = rhs.nRow;
    if (nCol*nRow > _DataSize) {
      free(_Data);
      _Data = (DataType *) malloc(nCol*nRow*sizeof(DataType));
      assert(_Data);
      _DataSize = nCol*nRow;
    }
    memcpy(_Data, rhs.getDataPtr(), nCol*nRow*sizeof(DataType));    
    _Transpose = rhs._Transpose;
    return *this;
  } 

  const Matrix<DataType>& operator= (const DataType *rhs) {
    if (sizeof(rhs)>_DataSize) {
      if (_Data) free(_Data);
      _Data = (DataType *) malloc( sizeof(rhs)*sizeof(DataType));
      assert(_Data);
      _DataSize = sizeof(rhs); 
    }
    memcpy(_Data, rhs, sizeof(rhs)*sizeof(DataType));
    _Transpose = CblasNoTrans;
    return *this;
  } 

  const Matrix<DataType>& operator= (const DataType rhs) {
    for (size_t i=0; i<_DataSize; i++) _Data[i]=rhs;
    return *this;
  } 

/* interface to private data */
  DataType * getDataPtr() const {
    return _Data;
  }

  size_t getDataSize() {
    return _DataSize;
  }
  void setTranspose(CBLAS_TRANSPOSE tr) {
    _Transpose = tr;
  }

  CBLAS_TRANSPOSE getTranspose() const {
    return _Transpose;
  }

  void setDim(size_t nr, size_t nc) {
    if (_DataSize < nc * nr) {
      DataType * newData = (DataType *) malloc( nc*nr*sizeof(DataType));
      assert(newData);
      memcpy(newData, _Data, nCol*nRow*sizeof(DataType));
      if (_Data) free(_Data);
      _Data = newData;
      _DataSize = nc*nr;
    }
    nCol = nc;
    nRow = nr;
  } 

/* operator overloading */
  inline DataType& operator()(size_t i, size_t j) {
    assert(i<nRow && j<nCol);
    if (_Transpose == CblasNoTrans) return _Data[j*nRow + i];
    else return _Data[i*nCol + j];
  }

  inline DataType operator()(size_t i, size_t j) const {
    assert(i<nRow && j<nCol);
    if (_Transpose == CblasNoTrans) return _Data[j*nRow + i];
    else return _Data[i*nCol + j];
  }

  // arithmetic
  const Matrix<DataType> operator*(const Matrix<DataType>& rhs) const {
    assert (nCol == rhs.nRow);
    size_t m,n,k1,k2,lda,ldb;
    m = nRow;
    k1 = nCol;
    k2 = rhs.nRow;
    n = rhs.nCol;
    lda = (_Transpose==CblasTrans)? nCol: nRow;
    ldb = (rhs._Transpose==CblasTrans)? rhs.nCol: rhs.nRow;

    Matrix<DataType> ma;

    ma.setDim(m,n);
    DataType alpha, beta;
    alpha = 1;
    beta = 0;
    
    gemm<DataType>(CblasColMajor, _Transpose, rhs._Transpose, m, n, k1, 
      alpha, _Data, lda, rhs.getDataPtr(), ldb, beta, ma.getDataPtr(), m);

    return ma;
  } 

  template<typename T> const Matrix<DataType> operator+(const Matrix<T>& rhs) {
    assert(nCol == rhs.nCol && nRow == rhs.nRow);
    Matrix<DataType> m(nRow, nCol);
    for (size_t i=0; i<nRow; i++)
      for (size_t j=0; j<nCol; j++)
        m(i,j) = (*this)(i,j)+rhs(i,j); 
    return m;
  }

  template<typename T> Matrix<DataType>& operator+=(const Matrix<T>& rhs) {
    assert(nCol == rhs.nCol && nRow == rhs.nRow);
    for (size_t i=0; i<nRow; i++)
      for (size_t j=0; j<nCol; j++)
        (*this)(i,j) += rhs(i,j); 
    return *this;
  }

  template<typename T> const Matrix<DataType> operator-(const Matrix<T>& rhs) {
    assert(nCol == rhs.nCol && nRow == rhs.nRow);
    Matrix<DataType> m(nRow, nCol);
    for (size_t i=0; i<nRow; i++)
      for (size_t j=0; j<nCol; j++)
        m(i,j) = (*this)(i,j)-rhs(i,j); 
    return m;
  }

  template<typename T> Matrix<DataType>& operator-=(const Matrix<T>& rhs) {
    assert(nCol == rhs.nCol && nRow == rhs.nRow);
    for (size_t i=0; i<nRow; i++)
      for (size_t j=0; j<nCol; j++)
        (*this)(i,j) -= rhs(i,j); 
    return *this;
  }

  template<typename T> const Matrix<DataType> operator/(const T rhs) {
    Matrix<DataType> m(nRow, nCol);
    for (size_t i=0; i<nRow; i++)
      for (size_t j=0; j<nCol; j++)
        m(i,j) = (*this)(i,j)/rhs; 
    return m;
  }

  template<typename T> Matrix<DataType>& operator/=(const T rhs) {
    for (size_t i=0; i<nRow; i++)
      for (size_t j=0; j<nCol; j++)
        (*this)(i,j) /= rhs; 
    return (*this);
  }

/* matrix specific operations */
  Matrix<DataType> &transpose() {
    _Transpose = (_Transpose == CblasTrans)? CblasNoTrans: CblasTrans;
   setDim(nCol, nRow);
   return *this;
  } 

  Matrix<DataType> &conjTrans() {
   transpose();
   for (size_t i=0; i<nRow; i++) 
     for (size_t j=0; j<nCol; j++) 
       (*this)(i,j)=conj((*this)(i,j));
   return *this;
  } 

protected:
  DataType * _Data; 
  size_t _DataSize;
  CBLAS_TRANSPOSE _Transpose;
};

/* Matrix derived classes */
template<typename DataType> class DMatrix : public Matrix<DataType> {
// diagnal matrix
public:
  DMatrix(size_t n) {
    Matrix<DataType>::_Data = (DataType *) malloc( n*n*sizeof(DataType));
    assert(Matrix<DataType>::_Data);
    Matrix<DataType>::nCol = n;
    Matrix<DataType>::nRow = n;
    Matrix<DataType>::_DataSize = n*n;
    
    for (size_t i=0; i<n; i++) 
      for (size_t j=0; j<n; j++)
        if (i==j) (*this)(i,j)=1;
        else (*this)(i,j)=0;
  }
};

template<> class DMatrix<complex<double> > : public Matrix<complex<double> > {
// complex diagnal
public:
  DMatrix(size_t n) {
    Matrix<complex<double> >::_Data = (complex<double>  *) malloc( n*n*sizeof(complex<double> ));
    assert(Matrix<complex<double> >::_Data);
    Matrix<complex<double> >::nCol = n;
    Matrix<complex<double> >::nRow = n;
    Matrix<complex<double> >::_DataSize = n*n;
    
    for (size_t i=0; i<n; i++) 
	    for (size_t j=0; j<n; j++)
        if (i==j) (*this)(i,j)=C1;
        else (*this)(i,j)=0;
  }
};


/* matrix operation functions */

template<typename DataType> Matrix<DataType> transpose(const Matrix<DataType> &m) {
  Matrix<DataType> m1 = m;
  m.transpose();
  return m1;
}

template<typename DataType> Matrix<DataType> conjTrans(const Matrix<DataType> &m) {
  Matrix<DataType> m1 = m;
  m.conjTrans();
  return m1;
}

/* arithmetic */
template<typename DataType, typename T> const Matrix<DataType> operator*(const T lhs, const Matrix<DataType> &ma) {
  Matrix<DataType> m(ma.nRow, ma.nCol);
  for (size_t i=0; i<ma.nRow; i++)
    for (size_t j=0; j<ma.nCol; j++)
      m(i,j) = ma(i,j)*lhs; 
  return m;
}

template<typename DataType> const Matrix<DataType> operator*(const Matrix<DataType>& lhs, const Matrix<DataType>& ma) {
  return ma.operator*(lhs);
}

/* stream operator overload */

template<typename DataType> ostream& operator<< (ostream& out, const Matrix<DataType> &m) {
  for (size_t i=0; i<m.nRow; i++) {
    for (size_t j=0; j<m.nCol; j++)
      out << m(i,j) << " ";
    if (i!=m.nRow-1) out << endl;
  }
  return out; 
}
 
}
#endif
