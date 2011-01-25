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

complex<double> C1 = complex<double>(1,0);
complex<double> Ci = complex<double>(0,1);

template<typename T> void gemm(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N, const MKL_INT K, const T alpha, const T *A, const MKL_INT lda, const T *B, const MKL_INT ldb, const T beta, T *c, const MKL_INT ldc) { };

template<> void gemm<double>(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N, const MKL_INT K, const double alpha, const double *A, const MKL_INT lda, const double *B, const MKL_INT ldb, const double beta, double *c, const MKL_INT ldc) { 
  cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, c, ldc); 
};

template<> void gemm<complex<double> >(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N, const MKL_INT K, const complex<double>  alpha, const complex<double>  *A, const MKL_INT lda, const complex<double>  *B, const MKL_INT ldb, const complex<double>  beta, complex<double>  *c, const MKL_INT ldc) { 
  cblas_zgemm(Order, TransA, TransB, M, N, K, &alpha, A, lda, B, ldb, &beta, c, ldc); 
};



template<typename DataType> class Matrix {
public:
  size_t nCol;
  size_t nRow;

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

  DataType * getDataPtr() const {
    return _Data;
  }

  size_t getDataSize() {
    return _DataSize;
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

  DataType& operator()(size_t i, size_t j) {
    assert(i<nRow && j<nCol);
    if (_Transpose == CblasNoTrans) return _Data[j*nRow + i];
    else return _Data[i*nCol + j];
  }

  DataType operator()(size_t i, size_t j) const {
    assert(i<nRow && j<nCol);
    if (_Transpose == CblasNoTrans) return _Data[j*nRow + i];
    else return _Data[i*nCol + j];
  }

  const Matrix<DataType> operator*(const Matrix<DataType>& rhs) {

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


template<typename DataType> class DMatrix : public Matrix<DataType> {

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

template<typename DataType> ostream& operator<< (ostream& out, const Matrix<DataType> &m) {
  for (size_t i=0; i<m.nRow; i++) {
    for (size_t j=0; j<m.nCol; j++)
      out << m(i,j) << " ";
    if (i!=m.nRow-1) out << endl;
  }
  return out; 
}
 
}

using namespace MAT;

int main() {

  Matrix<double> m1;

  double a[] = {1., 2., 3., 4.};
  double b[] = {1., 2., 3., 4., 5., 6.};

  m1 = a;
  m1.setDim(2, 2);
  std::cout << m1(0,1) << std::endl;

  Matrix<double> m2 = m1;
  m2(0,1) = 99.;

  Matrix<double> m3(b, 2, 3);

  cout << "m1=" << endl << m1 << endl;
  cout << "m2=" << endl << m2 << endl;
  cout << "m1*m2=" << endl <<m1*m2 << endl;

  cout << "m2=" << endl <<m3 << endl;
  cout << "m1*m1=" << endl <<m1*m1 << endl;

  Matrix<double> m4 = DMatrix<double>(3);

  cout << "m4=" << endl << m4 << endl;

  Matrix<double> m5(3,3);
  m5 = 3;

  cout << "m4*m5=" << endl << m4*m5 << endl;

  CMatrix m6(4,4);
  m6 = C1*3. + Ci*2.;

  cout << "m6=" << endl << m6 << endl;
  cout << "1=" << endl << CDMatrix(4) << endl;
  cout << "m6*1=" << endl << m6*CDMatrix(4) << endl;

  CD c1[]={ CD(1.,0.), CD(-3., 0.), CD(2.0, 0.), CD(-2., .0), CD(6., 0.), CD(3.0, 0.)};
  CMatrix m7(c1, 2,3);
  //cout << "t(m7)=" << endl << transpose(m7) << endl;
  //cout << "t^2(m7)=" << endl << transpose(transpose(m7)) << endl;
  cout << "m7=" << endl << m7 << endl;

  CD c2[]={CD( 1., 1.0),  CD( 2., 0.),    CD(3., 0.0 ),   CD(6., 0.), 
           CD( 1., .0),  CD( 2., 0.),    CD(3., 0.0 ),   CD(6., 0.), 
           CD( 1., .0),  CD( 2., 0.),    CD(3., 0.0 ),   CD(6., 0.)};
  CMatrix m8(c2, 4, 3);
  m8.conjTrans();
  cout << "m8=" << endl << m8 << endl;

  cout << "m7*m8=" << endl << m7*m8 << endl;
  return 0;
}
