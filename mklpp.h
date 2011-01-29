#ifndef MKLPP_H
#include <iostream>
#include <mkl_cblas.h>
#include <mkl_lapacke.h>
#include <mkl.h>
#include <string.h>
#include <assert.h>
#include <memory>

namespace mklpp {
using namespace std;

//#define CD MKL_Complex16 
#define CD ComplexDouble 
#define MKLCD MKL_Complex16 
#define cmatrix matrix<CD>

/* wrapping MKL_Complex16 */
struct _ComplexDouble : MKLCD {
  _ComplexDouble() {};
  _ComplexDouble(double r, double i) { real = r; imag = i;};
  _ComplexDouble(double r) { real = r; imag =0.0;};
  _ComplexDouble(MKLCD c) { real = c.real; imag =c.imag;};

  const struct _ComplexDouble & operator=(double r) {
    real = r;
    imag = 0.0;
    return *this;
  }

  const struct _ComplexDouble & operator=(const MKLCD &c) {
    real = c.real;
    imag = c.imag;
    return *this;
  }
};

typedef struct _ComplexDouble ComplexDouble;

/* complex constant */
CD C1 (1.,0.);
CD Ci (0.,1.);
CD C0 (0.,0.);

/* overloading MKL_Complex16 operations */
inline const CD operator+(const CD &a, const CD &b) {
  CD c;
  c.real = a.real + b.real; 
  c.imag = a.imag + b.real;
  return c;
}

inline CD& operator+=(CD &a, const CD &b) {
  a.real += b.real; 
  a.imag += b.real;
  return a;
}

inline const CD operator-(const CD &a, const CD &b) {
  CD c;
  c.real = a.real - b.real; 
  c.imag = a.imag - b.imag;
  return c;
}

inline CD& operator-=(CD &a, const CD &b) {
  a.real -= b.real; 
  a.imag -= b.real;
  return a;
}

inline const CD operator*(const CD &a, const CD &b) {
  CD c;
  c.real = a.real*b.real - a.imag*b.imag; 
  c.imag = a.real*b.imag + a.imag*b.real;
  return c;
}

inline CD& operator*=(CD &a, const CD &b) {
  a.real = a.real*b.real - a.imag*b.imag; 
  a.imag = a.real*b.imag + a.imag*b.real;
  return a;
}

inline const CD operator*(const CD &a, const double b) {
  CD c;
  c.real = a.real * b; 
  c.imag = a.imag * b;
  return c;
}

inline CD& operator*=(CD &a, const double b) {
  a.real = a.real * b; 
  a.imag = a.imag * b;
  return a;
}

inline const CD operator*(const double b, const CD &a) {
  CD c;
  c.real = a.real * b; 
  c.imag = a.imag * b;
  return c;
}

inline const CD operator/(const CD &a, const double b) {
  CD c;
  c.real = a.real / b; 
  c.imag = a.imag / b;
  return c;
}

inline const CD conj(const CD &a) {
  CD c;
  c.real = a.real; 
  c.imag = -a.imag;
  return c;
}

ostream& operator<< (ostream& out, const CD &m) {
  out << "(" << m.real << ", " << m.imag << ")";
  return out; 
}
 

/* MKL function wrappers implemented with templates*/
template<typename T> void gemm(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N, const MKL_INT K, const T *A, const MKL_INT lda, const T *B, const MKL_INT ldb, T *c, const MKL_INT ldc) { 
  assert(0); // always fails
};

template<> void gemm<double>(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N, const MKL_INT K, const double *A, const MKL_INT lda, const double *B, const MKL_INT ldb, double *c, const MKL_INT ldc) { 
  cblas_dgemm(Order, TransA, TransB, M, N, K, 1., A, lda, B, ldb, 0., c, ldc); 
};

template<> void gemm<CD>(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N, const MKL_INT K, const CD  *A, const MKL_INT lda, const CD  *B, const MKL_INT ldb, CD  *c, const MKL_INT ldc) { 
  CD alpha (1., 0.);
  CD beta (0., 0.);
  cblas_zgemm(Order, TransA, TransB, M, N, K, &alpha, A, lda, B, ldb, &beta, c, ldc); 
};

template<typename T> int geev(const int order, char nl, char nr, const MKL_INT n, T *a, const MKL_INT lda, T * w, T *vl, const MKL_INT ldvl, T *vr, const MKL_INT ldvr) { 
  assert(0); // always fails
};

template<> int geev<CD>(const int order, char nl, char nr, const MKL_INT n, CD *a, const MKL_INT lda, CD * w, CD *vl, const MKL_INT ldvl, CD *vr, const MKL_INT ldvr) { 
  return LAPACKE_zgeev(LAPACK_COL_MAJOR, nl, nr, n, a, lda, w, vl, ldvl, vr, ldvr); 
};

/* array for auto_ptr */
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

  template<typename T> DataArray(const T *da, const size_t s) {
    _Data = (DataType *) malloc(s*sizeof(DataType));
    assert (_Data);
    if (da) for (size_t i=0; i<s; i++) _Data[i] = da[i];
    size = s;
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

  template<typename T> matrix(const T *data, size_t nr, size_t nc) : nCol(nc), nRow(nr) {
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

  inline DataType& operator()(size_t i) {
    assert((nCol == 1 && i<nRow) || (nRow == 1 && i<nCol));
    return (*_Data)[i];
  }

  inline DataType& operator()(size_t i) const {
    assert((nCol == 1 && i<nRow) || (nRow == 1 && i<nCol));
    return (*_Data)[i];
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

    gemm<DataType>(CblasColMajor, _Transpose, rhs._Transpose, m, n, k1, 
      getDataPtr(), lda, rhs.getDataPtr(), ldb, ma.getDataPtr(), m);

    return ma;
  } 

  template<typename T> const matrix<DataType> operator+(const matrix<T>& rhs) const {
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

  template<typename T> const matrix<DataType> operator-(const matrix<T>& rhs) const {
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
       (*this)(i,j)=mklpp::conj((*this)(i,j));
   return *this;
  } 

  matrix<DataType> &conj() {
   for (size_t i=0; i<nRow; i++) 
     for (size_t j=0; j<nCol; j++) 
       (*this)(i,j)=mklpp::conj((*this)(i,j));
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

/* solvers */
  void eigen(matrix<DataType> &e, matrix<DataType> &vl, matrix<DataType> &vr) {
    assert (nCol == nRow && (_Transpose == CblasNoTrans || _Transpose == CblasTrans));
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

    geev<DataType>(LAPACK_COL_MAJOR, nl, nr, nCol, m.getDataPtr(), nCol, e.getDataPtr(), vl.getDataPtr(), nCol, vr.getDataPtr(), nCol); 
  } 

  void reigen(matrix<DataType> &e, matrix<DataType> &vr) {
    assert (nCol == nRow && (_Transpose == CblasNoTrans || _Transpose == CblasTrans));
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

    geev<DataType>(LAPACK_COL_MAJOR, nl, nr, nCol, m.getDataPtr(), nCol, e.getDataPtr(), NULL, nCol, vr.getDataPtr(), nCol); 
  } 

  void leigen(matrix<DataType> &e, matrix<DataType> &vl) {
    assert (nCol == nRow && (_Transpose == CblasNoTrans || _Transpose == CblasTrans));
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

    geev<DataType>(LAPACK_COL_MAJOR, nl, nr, nCol, m.getDataPtr(), nCol, e.getDataPtr(), vl.getDataPtr(), nCol, NULL, nCol); 
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
    matrix<DataType>::_Transpose = CblasNoTrans;
    
    for (size_t i=0; i<n; i++) 
      for (size_t j=0; j<n; j++)
        if (i==j) (*this)(i,j)=1.;
        else (*this)(i,j)=0.;
  }
};

typedef idmatrix<CD> cidmatrix;

template<typename DataType> class diagmatrix : public matrix<DataType> {
// diagnol matrix build from vector
public:
  template<typename T> diagmatrix(const matrix<T> &v) {
    assert (v.nCol == 1 || v.nRow == 1); 
    size_t n = (v.nCol >v.nRow? v.nCol:v.nRow);
    matrix<DataType>::_Data = DataPtr<DataType>::Type(new DataArray<DataType>(n*n));
    matrix<DataType>::nCol = n;
    matrix<DataType>::nRow = n;
    matrix<DataType>::_DataSize = n*n;
    matrix<DataType>::_Transpose = CblasNoTrans;
    
    for (size_t i=0; i<n; i++) 
      for (size_t j=0; j<n; j++)
        if (i==j) (*this)(i,j)=v(i);
        else (*this)(i,j)=0.;
  }
};

typedef diagmatrix<CD> cdiagmatrix;

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

template<typename DataType, typename T> const matrix<DataType> & operator*=(matrix<DataType> &ma, const T lhs) {
  for (size_t i=0; i<ma.nRow; i++)
    for (size_t j=0; j<ma.nCol; j++)
      ma(i,j) *= lhs; 
  return ma;
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
