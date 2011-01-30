#ifndef MKLPP_MATRIX_H
#define MKLPP_MATRIX_H
#include "mklpp_complex.h"
#include "mklpp_wrapper.h"
#include "mklpp_dataarray.h"

namespace mklpp {

/////////////////////////////////////////////////////////////////////////////
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

  template<typename T> matrix(const matrix<T> &m) : nCol(m.nCol), nRow(m.nRow), _DataSize(m.nRow*m.nCol), _Transpose(m.getTranspose()){
    _Data = DataPtr<DataType>::Type(new DataArray<DataType>(m.getDataPtr(), m.nRow*m.nCol)); 
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

  DataArray<DataType>& getDataArray() const {
    return *_Data;
  }

  size_t getDataSize() const {
    return _DataSize;
  }

  CBLAS_TRANSPOSE getTranspose() const {
    return _Transpose;
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

  inline DataType& operator[](size_t i) {
    assert((nCol == 1 && i<nRow) || (nRow == 1 && i<nCol));
    return (*_Data)[i];
  }

  inline DataType& operator[](size_t i) const {
    assert((nCol == 1 && i<nRow) || (nRow == 1 && i<nCol));
    return (*_Data)[i];
  }

  // arithmetic

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

  const matrix<DataType> operator*(const double rhs) {
    matrix<DataType> m(nRow, nCol);
    for (size_t i=0; i<nRow; i++)
      for (size_t j=0; j<nCol; j++)
        m(i,j) = (*this)(i,j)*rhs; 
    return m;
  }

  const matrix<DataType> operator*(const CD rhs) {
    matrix<DataType> m(nRow, nCol);
    for (size_t i=0; i<nRow; i++)
      for (size_t j=0; j<nCol; j++)
        m(i,j) = (*this)(i,j)*rhs; 
    return m;
  }

  matrix<DataType> & operator*=(const double lhs) {
    for (size_t i=0; i<ma.nRow; i++)
      for (size_t j=0; j<ma.nCol; j++)
      (*this)(i,j) *= lhs; 
    return *this;
  }

  matrix<DataType> & operator*=(const CD lhs) {
    for (size_t i=0; i<ma.nRow; i++)
      for (size_t j=0; j<ma.nCol; j++)
      (*this)(i,j) *= lhs; 
    return *this;
  }
 
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

  matrix<DataType>& operator*=(const matrix<DataType>& rhs) {
    matrix<DataType> ma = (*this)*rhs; 
    (*this) = ma;
    return (*this);
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

typedef matrix<CD> cmatrix;
typedef matrix<double> dmatrix;

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

template<typename DataType> matrix<DataType> trans(const matrix<DataType> &m) {
  matrix<DataType> m1 = m;
  m.trans();
  return m1;
}

template<typename DataType> matrix<DataType> herm(const matrix<DataType> &m) {
  matrix<DataType> m1 = m;
  m1.herm();
  return m1;
}

template<typename DataType> matrix<DataType> conj(const matrix<DataType> &m) {
  matrix<DataType> m1 = m;
  m.conj();
  return m1;
}

/* arithmetic */
template<typename DataType> const matrix<DataType> operator*(const double lhs, const matrix<DataType> &ma) {
  matrix<DataType> m(ma.nRow, ma.nCol);
  for (size_t i=0; i<ma.nRow; i++)
    for (size_t j=0; j<ma.nCol; j++)
      m(i,j) = ma(i,j)*lhs; 
  return m;
}

template<typename DataType> const matrix<DataType> operator*(const CD lhs, const matrix<DataType> &ma) {
  matrix<DataType> m(ma.nRow, ma.nCol);
  for (size_t i=0; i<ma.nRow; i++)
    for (size_t j=0; j<ma.nCol; j++)
      m(i,j) = ma(i,j)*lhs; 
  return m;
}

/* stream operator overload */
template<typename DataType> std::ostream& operator<< (std::ostream& out, const matrix<DataType> &m) {
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
template<typename DataType> class idmatrix : public matrix<DataType> {
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
typedef idmatrix<double> didmatrix;
/////////////////////////////////////////////////////////////////////////////


}
#endif
