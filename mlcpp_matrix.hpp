#ifndef MKLPP_MATRIX_HPP
#define MKLPP_MATRIX_HPP
#include <mlcpp_complex.hpp>
#include <mlcpp_wrapper.hpp>
#include <mlcpp_dataarray.hpp>
#include <sstream>

namespace mlcpp {

/////////////////////////////////////////////////////////////////////////////
/* matrix class */
template<typename DataType> 
class matrix {
public:
  ~matrix() {}

  matrix();

  matrix(size_t nr, size_t nc);
 
  template<typename T> 
  matrix(const T *data, size_t nr, size_t nc); 

  template<typename T> 
  matrix(const matrix<T> &m);
 
  matrix(const matrix<DataType> &m); 

  matrix<DataType>& operator= (const matrix<DataType> &rhs);
 
  matrix<DataType>& operator= (const DataType rhs);
 
  bool operator== (const matrix<DataType> &rhs); 

  bool operator!= (const matrix<DataType> &rhs); 

/* interface to private data */
  size_t nRow() const;

  size_t nCol() const;

  DataType * getDataPtr() const; 

  DataArray<DataType>& getDataArray() const; 

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

  const matrix<DataType> operator*(const double rhs) const;

  const matrix<DataType> operator*(const CD rhs) const;

  matrix<DataType> & operator*=(const double rhs);

  matrix<DataType> & operator*=(const CD rhs);
 
  const matrix<DataType> operator*(const matrix<DataType>& rhs) const;

  matrix<DataType>& operator*=(const matrix<DataType>& rhs);
 
  template<typename T> const matrix<DataType> operator/(const T rhs) const;

  template<typename T> matrix<DataType>& operator/=(const T rhs);

/* matrix specific operations */
  matrix<DataType> &trans();

  matrix<DataType> &herm();

  matrix<DataType> &conj();

  matrix<DataType> block(size_t r1, size_t r2, size_t c1, size_t c2) const;

  matrix<DataType> &replace(size_t r, size_t c, const matrix<DataType> &mt);

  matrix<DataType> col(size_t c) const;

  matrix<DataType> row(size_t r) const;

/* solvers */
  void eigen(matrix<CD> &e, matrix<DataType> &vl, matrix<DataType> &vr);

  void reigen(matrix<CD> &e, matrix<DataType> &vr);

  void leigen(matrix<CD> &e, matrix<DataType> &vl);

/* private data */
protected:
  size_t m_ncol;

  size_t m_nrow;

  typename DataPtr<DataType>::Type m_data; 

  CBLAS_TRANSPOSE m_transpose;
};

typedef matrix<CD> cmatrix;
typedef matrix<double> dmatrix;

template<typename DataType> 
matrix<DataType>::matrix() : m_ncol(0), m_nrow(0) {
  m_transpose = CblasNoTrans;
} 

template<typename DataType> 
matrix<DataType>::matrix(size_t nr, size_t nc) : m_ncol(nc), m_nrow(nr) {
    m_data = typename DataPtr<DataType>::Type(new DataArray<DataType>(m_ncol*m_nrow));
    m_transpose = CblasNoTrans;
  } 

template<typename DataType> 
template<typename T> 
matrix<DataType>::matrix(const T *data, size_t nr, size_t nc) 
  : m_ncol(nc), m_nrow(nr) {
  m_data =typename DataPtr<DataType>::Type(new DataArray<DataType>(data, nr*nc));
  m_transpose = CblasNoTrans;
} 

template<typename DataType> 
template<typename T> 
matrix<DataType>::matrix(const matrix<T> &m) 
  : m_ncol(m.nCol()), m_nrow(m.nRow()), 
  m_transpose(m.getTranspose()){
  m_data =typename DataPtr<DataType>::Type(
    new DataArray<DataType>(m.getDataPtr(), m.nRow()*m.nCol())
  ); 
  m_transpose = m.getTranspose();
} 

template<typename DataType> 
matrix<DataType>::matrix(const matrix<DataType> &m) 
  : m_ncol(m.nCol()), m_nrow(m.nRow()), 
  m_transpose(m.m_transpose){
  m_data =typename DataPtr<DataType>::Type(new DataArray<DataType>(*(m.m_data)));
}

template<typename DataType> 
matrix<DataType>& 
matrix<DataType>::operator= (const matrix<DataType> &rhs) {
  if (&rhs==this) return *this;
  m_ncol = rhs.nCol();
  m_nrow = rhs.nRow();
  m_data =typename DataPtr<DataType>::Type(new DataArray<DataType>(*(rhs.m_data)));
  m_transpose = rhs.m_transpose;
  return *this;
} 

template<typename DataType> 
matrix<DataType>& 
matrix<DataType>::operator= (const DataType rhs) {
  for (size_t i=0; i<m_ncol*m_nrow; i++) (*m_data)[i]=rhs;
  return *this;
} 

template<typename DataType> 
bool matrix<DataType>::operator== (const matrix<DataType> &rhs) {
  if (&rhs==this) return true;
  if (m_nrow != rhs.nRow() || m_ncol != rhs.nCol()) return false;
  for (size_t i=0; i<m_nrow; i++)
    for (size_t j=0; j<m_ncol; j++)
      if ((*this)(i,j)!=rhs(i,j)) return false;
  return true;
} 

template<typename DataType> 
bool matrix<DataType>::operator!= (const matrix<DataType> &rhs) {
  return !(operator==(rhs));
}

/* interface to private data */
template<typename DataType> 
size_t matrix<DataType>::nRow() const {
  return m_nrow;
}

template<typename DataType> 
size_t matrix<DataType>::nCol() const {
  return m_ncol;
}

template<typename DataType> 
DataType * matrix<DataType>::getDataPtr() const {
  return (*m_data).getDataPtr();
}

template<typename DataType> 
DataArray<DataType>& matrix<DataType>::getDataArray() const {
  return *m_data;
}

template<typename DataType> 
CBLAS_TRANSPOSE matrix<DataType>::getTranspose() const {
  return m_transpose;
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

/* operator overloading */
template<typename DataType> 
inline DataType& matrix<DataType>::operator()(size_t i, size_t j) {
  assert(i<m_nrow && j<m_ncol);
  if (m_transpose == CblasNoTrans) return (*m_data)[j*m_nrow + i];
  else return (*m_data)[i*m_ncol + j];
}

template<typename DataType> 
inline DataType matrix<DataType>::operator()(size_t i, size_t j) const {
  assert(i<m_nrow && j<m_ncol);
  if (m_transpose == CblasNoTrans) return (*m_data)[j*m_nrow + i];
  else return (*m_data)[i*m_ncol + j];
}

template<typename DataType> 
inline DataType& matrix<DataType>::operator()(size_t i) {
  return (*m_data)[i];
}

template<typename DataType> 
inline DataType& matrix<DataType>::operator()(size_t i) const {
  return (*m_data)[i];
}

template<typename DataType> 
inline DataType& matrix<DataType>::operator[](size_t i) {
  return (*m_data)[i];
}

template<typename DataType> 
inline DataType& matrix<DataType>::operator[](size_t i) const {
  return (*m_data)[i];
}

// arithmetic

template<typename DataType> 
template<typename T> 
const matrix<DataType> 
matrix<DataType>::operator+(const matrix<T>& rhs) const {
  assert(m_ncol == rhs.nCol() && m_nrow == rhs.nRow());
  matrix<DataType> m(m_nrow, m_ncol);
  for (size_t i=0; i<m_nrow; i++)
    for (size_t j=0; j<m_ncol; j++)
      m(i,j) = (*this)(i,j)+rhs(i,j); 
  return m;
}

template<typename DataType> 
template<typename T> 
matrix<DataType>& matrix<DataType>::operator+=(const matrix<T>& rhs) {
  assert(m_ncol == rhs.nCol() && m_nrow == rhs.nRow());
  for (size_t i=0; i<m_nrow; i++)
    for (size_t j=0; j<m_ncol; j++)
      (*this)(i,j) += rhs(i,j); 
  return *this;
}

template<typename DataType> 
template<typename T>
const matrix<DataType> 
matrix<DataType>::operator-(const matrix<T>& rhs) const {
  assert(m_ncol == rhs.nCol() && m_nrow == rhs.nRow());
  matrix<DataType> m(m_nrow, m_ncol);
  for (size_t i=0; i<m_nrow; i++)
    for (size_t j=0; j<m_ncol; j++)
      m(i,j) = (*this)(i,j)-rhs(i,j); 
  return m;
}

template<typename DataType> 
template<typename T> 
matrix<DataType>& matrix<DataType>::operator-=(const matrix<T>& rhs) {
  assert(m_ncol == rhs.nCol() && m_nrow == rhs.nRow());
  for (size_t i=0; i<m_nrow; i++)
    for (size_t j=0; j<m_ncol; j++)
      (*this)(i,j) -= rhs(i,j); 
  return *this;
}

template<typename DataType> 
const matrix<DataType> matrix<DataType>::operator*(const double rhs) const {
  matrix<DataType> m(m_nrow, m_ncol);
  for (size_t i=0; i<m_nrow; i++)
    for (size_t j=0; j<m_ncol; j++)
      m(i,j) = (*this)(i,j)*rhs; 
  return m;
}

template<typename DataType> 
const matrix<DataType> matrix<DataType>::operator*(const CD rhs) const {
  matrix<DataType> m(m_nrow, m_ncol);
  for (size_t i=0; i<m_nrow; i++)
    for (size_t j=0; j<m_ncol; j++)
      m(i,j) = (*this)(i,j)*rhs; 
  return m;
}

template<typename DataType> 
matrix<DataType> & matrix<DataType>::operator*=(const double rhs) {
  for (size_t i=0; i<nRow(); i++)
    for (size_t j=0; j<nCol(); j++)
    (*this)(i,j) *= rhs; 
  return *this;
}

template<typename DataType> 
matrix<DataType> & matrix<DataType>::operator*=(const CD rhs) {
  for (size_t i=0; i<nRow(); i++)
    for (size_t j=0; j<nCol(); j++)
    (*this)(i,j) *= rhs; 
  return *this;
}
 
template<typename DataType> 
const matrix<DataType> 
matrix<DataType>::operator*(const matrix<DataType>& rhs) const {
  assert (m_ncol == rhs.nRow());
  size_t m,n,k1,k2,lda,ldb;
  m = m_nrow;
  k1 = m_ncol;
  k2 = rhs.nRow();
  n = rhs.nCol();
  lda = (m_transpose==CblasTrans)? m_ncol: m_nrow;
  ldb = (rhs.m_transpose==CblasTrans)? rhs.nCol(): rhs.nRow();

  matrix<DataType> ma;

  ma.resize(m,n);

  gemm<DataType>(CblasColMajor, m_transpose, rhs.m_transpose, m, n, k1, 
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
  for (size_t i=0; i<m_nrow; i++)
    for (size_t j=0; j<m_ncol; j++)
      m(i,j) = (*this)(i,j)/rhs; 
  return m;
}

template<typename DataType> 
template<typename T> 
matrix<DataType>& matrix<DataType>::operator/=(const T rhs) {
  for (size_t i=0; i<m_nrow; i++)
    for (size_t j=0; j<m_ncol; j++)
      (*this)(i,j) /= rhs; 
  return (*this);
}

/* matrix specific operations */
template<typename DataType> 
matrix<DataType> & matrix<DataType>::trans() {
  m_transpose = (m_transpose == CblasTrans)? CblasNoTrans: CblasTrans;
 resize(m_ncol, m_nrow);
 return *this;
} 

template<typename DataType> 
matrix<DataType> & matrix<DataType>::herm() {
 trans();
 for (size_t i=0; i<m_nrow; i++) 
   for (size_t j=0; j<m_ncol; j++) 
     (*this)(i,j)=mlcpp::conj((*this)(i,j));
 return *this;
} 

template<typename DataType> 
matrix<DataType> & matrix<DataType>::conj() {
 for (size_t i=0; i<m_nrow; i++) 
   for (size_t j=0; j<m_ncol; j++) 
     (*this)(i,j)=mlcpp::conj((*this)(i,j));
 return *this;
} 

template<typename DataType> 
matrix<DataType> matrix<DataType>::block(size_t r1, size_t r2, 
                                         size_t c1, size_t c2) const {
  assert(r1<=r2 && c1<=c2);
  matrix<DataType> m(r2-r1, c2-c1);
  for (size_t i=0; i<m.nRow(); i++) 
    for (size_t j=0; j<m.nCol(); j++) 
      m(i,j) = (*this)(i+r1, j+c1);
  return m;
}

template<typename DataType> 
matrix<DataType> & matrix<DataType>::replace(size_t r, size_t c, 
                                            const matrix<DataType> &mt) {
  for (size_t i=0; i<mt.nRow(); i++) 
    for (size_t j=0; j<mt.nCol(); j++) 
      (*this)(i+r, j+c) = mt(i,j);
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

/* solvers */
template<typename DataType> 
void matrix<DataType>::eigen(matrix<CD> &e, matrix<DataType> &vl, 
                             matrix<DataType> &vr) {
  assert (m_ncol == m_nrow && 
          (m_transpose == CblasNoTrans || m_transpose == CblasTrans));
  matrix<DataType> m(*this);
  if ( m_transpose == CblasTrans) {
    m.resize(m_ncol, m_ncol);
    for (size_t i=0; i<m_nrow; i++) 
      for (size_t j=0; j<m_ncol; j++) 
        m(i,j)= (*this)(j,i);
  }    

  vl.resize(m_ncol,m_ncol);
  vr.resize(m_ncol,m_ncol);
  e.resize(m_ncol, 1);
  
  char nl = &vl? 'V': 'N';
  char nr = &vr? 'V': 'N';

  geev<DataType>(LAPACK_COL_MAJOR, nl, nr, m_ncol, m.getDataPtr(), m_ncol, 
                 e.getDataPtr(), vl.getDataPtr(), m_ncol, 
                 vr.getDataPtr(), m_ncol);
} 

template<typename DataType> 
void matrix<DataType>::reigen(matrix<CD> &e, matrix<DataType> &vr) {
  assert (m_ncol == m_nrow && (m_transpose == CblasNoTrans || 
                           m_transpose == CblasTrans));
  matrix<DataType> m(*this);
  if ( m_transpose == CblasTrans) {
    m.resize(m_ncol, m_ncol);
    for (size_t i=0; i<m_nrow; i++) 
      for (size_t j=0; j<m_ncol; j++) 
        m(i,j) = (*this)(j,i);
  }    

  vr.resize(m_ncol,m_ncol);
  e.resize(m_ncol, 1);
  
  char nl = 'N';
  char nr = 'V';

  geev<DataType>(LAPACK_COL_MAJOR, nl, nr, m_ncol, m.getDataPtr(), m_ncol, 
                 e.getDataPtr(), NULL, m_ncol, vr.getDataPtr(), m_ncol); 
} 

template<typename DataType> 
void matrix<DataType>::leigen(matrix<CD> &e, matrix<DataType> &vl) {
  assert (m_ncol == m_nrow && (m_transpose == CblasNoTrans || 
                           m_transpose == CblasTrans));
  matrix<DataType> m(*this);
  if ( m_transpose == CblasTrans) {
    m.resize(m_ncol, m_ncol);
    for (size_t i=0; i<m_nrow; i++) 
      for (size_t j=0; j<m_ncol; j++) 
        m(i,j) = (*this)(j,i);
  }    

  vl.resize(m_ncol,m_ncol);
  e.resize(m_ncol, 1);
  
  char nl = 'V';
  char nr = 'N';

  geev<DataType>(LAPACK_COL_MAJOR, nl, nr, m_ncol, m.getDataPtr(), m_ncol, 
                 e.getDataPtr(), vl.getDataPtr(), m_ncol, NULL, m_ncol); 
} 
/////////////////////////////////////////////////////////////////////////////

matrix<double> real(const matrix<CD> &ma) {
 matrix<double> m(ma.nRow(), ma.nCol());
 for (size_t i=0; i<ma.nRow(); i++) 
   for (size_t j=0; j<ma.nCol(); j++) 
     m(i,j)=ma(i,j).real;
 return m;
} 

matrix<double> imag(const matrix<CD> &ma) {
 matrix<double> m(ma.nRow(), ma.nCol());
 for (size_t i=0; i<ma.nRow(); i++) 
   for (size_t j=0; j<ma.nCol(); j++) 
     m(i,j)=ma(i,j).imag;
 return m;
}

template<typename DataType> 
matrix<DataType> trans(const matrix<DataType> &m) {
  matrix<DataType> m1 = m;
  m1.trans();
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
  m1.conj();
  return m1;
}

/* arithmetic */
template<typename DataType> 
const matrix<DataType> operator*(const double lhs, 
                                 const matrix<DataType> &ma) {
  matrix<DataType> m(ma.nRow(), ma.nCol());
  for (size_t i=0; i<ma.nRow(); i++)
    for (size_t j=0; j<ma.nCol(); j++)
      m(i,j) = ma(i,j)*lhs; 
  return m;
}

template<typename DataType> 
const matrix<DataType> operator*(const CD lhs, 
                                 const matrix<DataType> &ma) {
  matrix<DataType> m(ma.nRow(), ma.nCol());
  for (size_t i=0; i<ma.nRow(); i++)
    for (size_t j=0; j<ma.nCol(); j++)
      m(i,j) = ma(i,j)*lhs; 
  return m;
}

/* stream operator overload */
template<typename DataType> 
std::ostream& operator<< (std::ostream& out, 
                          const matrix<DataType> &m) {
  for (size_t i=0; i<m.nRow(); i++) {
    for (size_t j=0; j<m.nCol(); j++)
      out << m(i,j) << " ";
    if (i!=m.nRow()-1) out << std::endl;
  }
  return out; 
}

template<typename DataType> 
std::string toString(const matrix<DataType> &m) {
  std::ostringstream out;
  out << m;
  return out.str(); 
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
  matrix<DataType>::m_data = 
    typename DataPtr<DataType>::Type(new DataArray<DataType>(n*n));
  matrix<DataType>::m_ncol = n;
  matrix<DataType>::m_nrow = n;
  matrix<DataType>::m_transpose = CblasNoTrans;
  
  for (size_t i=0; i<n; i++) 
    for (size_t j=0; j<n; j++)
      if (i==j) (*this)(i,j)=1.;
      else (*this)(i,j)=0.;
}
/////////////////////////////////////////////////////////////////////////////


}
#endif
