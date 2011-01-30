#ifndef MKLPP_VECTOR_H
#define MKLPP_VECTOR_H

#include <mklpp_matrix.hpp>
#include <algorithm>
#include <vector>

namespace mklpp {

/////////////////////////////////////////////////////////////////////////////
/* vector class
   nCol = 1 or nRow =1 matrix */
template<typename DataType> class vector : public matrix<DataType> {
public:
  vector();

  vector(size_t i);

  vector(const matrix<DataType> &m);

  template<typename T> 
  vector(const T *data, size_t n);

  vector(const std::vector<DataType> &v);

  vector<DataType>& operator=(const vector<DataType> &v);

  vector<DataType>& operator=(const std::vector<DataType> &v);

  void resize(size_t n);
  void resize(size_t nc, size_t nr);

  size_t size() const;

  template<typename T> 
  const vector<DataType> operator+(const vector<T>& rhs) const;
 
  template<typename T> 
  vector<DataType>& operator+=(const vector<T>& rhs);
 
  template<typename T> 
  const vector<DataType> operator-(const vector<T>& rhs) const;
 
  template<typename T> 
  vector<DataType>& operator-=(const vector<T>& rhs); 

  const vector<DataType> operator*(const double rhs);

  vector<DataType>& operator*=(const double rhs); 

  const vector<DataType> operator*(const CD rhs);

  vector<DataType>& operator*=(const CD rhs);

  const vector<DataType> operator*(const matrix<DataType>& rhs) const;

  vector<DataType>& operator*=(const matrix<DataType>& rhs);

  template<typename T> 
  const DataType operator* (const vector<T> &rv);

  template<typename T> 
  const vector<DataType> operator/(const T rhs) const;

  template<typename T> 
  vector<DataType>& operator/=(const T rhs); 

  vector<DataType> &trans();

  vector<DataType> &herm();

  vector<DataType> &conj(); 

  vector<DataType> block(size_t i1, size_t i2);

  vector<DataType> &insert(size_t i, vector<DataType> v);

  vector<DataType> &sort();

  std::vector<DataType> stl();
};

typedef vector<CD> cvector;
typedef vector<double> dvector;

template<typename DataType>
vector<DataType>::vector() : matrix<DataType>() { }

template<typename DataType>
vector<DataType>::vector(size_t i) : matrix<DataType>(i, 1) { }

template<typename DataType>
vector<DataType>::vector(const matrix<DataType> &m) {
  matrix<DataType>::_Data = 
    DataPtr<DataType>::Type(new DataArray<DataType>(m.getDataPtr(), 
                                                    m.nRow*m.nCol)); 
  matrix<DataType>::nCol =1;
  matrix<DataType>::nRow = m.nRow*m.nCol;
  matrix<DataType>::_DataSize = m.nRow*m.nCol;
  matrix<DataType>::_Transpose = CblasNoTrans;
}

template<typename DataType>
template<typename T> 
vector<DataType>::vector(const T *data, size_t n) 
  : matrix<DataType>(data, n, 1){ } 

template<typename DataType>
vector<DataType>::vector(const std::vector<DataType> &v) 
  : matrix<DataType>(&v[0], v.size(), 1) {
}

template<typename DataType>
vector<DataType>& vector<DataType>::operator=(const vector<DataType> &v) { 
  matrix<DataType>::_Data = 
    DataPtr<DataType>::Type(new DataArray<DataType>(v.getDataPtr(), 
                                                    v.size())); 
  matrix<DataType>::nCol =1;
  matrix<DataType>::nRow = v.size();
  matrix<DataType>::_DataSize = v.size();
  matrix<DataType>::_Transpose = CblasNoTrans;
  return *this;
} 

template<typename DataType>
vector<DataType> & 
vector<DataType>::operator=(const std::vector<DataType> &v) { 
  resize(v.size());
  for(size_t i=0; i<v.size(); i++) (*this)[i] = v[i];
  return *this;
} 

template<typename DataType>
void vector<DataType>::resize(size_t n) { matrix<DataType>::resize(n,1);}

template<typename DataType>
void vector<DataType>::resize(size_t nc, size_t nr) {};

template<typename DataType>
size_t vector<DataType>::size() const {
  return matrix<DataType>::nCol==1?matrix<DataType>::nRow:matrix<DataType>::nCol;
}

template<typename DataType>
template<typename T> 
const vector<DataType> vector<DataType>::operator+(const vector<T>& rhs) const {
  return matrix<DataType>::operator+(rhs);
}

template<typename DataType>
template<typename T> 
vector<DataType>& vector<DataType>::operator+=(const vector<T>& rhs) {
  return matrix<DataType>::operator+=(rhs);
}

template<typename DataType>
template<typename T> 
const vector<DataType> 
vector<DataType>::operator-(const vector<T>& rhs) const {
  return matrix<DataType>::operator-(rhs);
}

template<typename DataType>
template<typename T> 
vector<DataType>& vector<DataType>::operator-=(const vector<T>& rhs) {
  return matrix<DataType>::operator-=(rhs);
}

template<typename DataType>
const vector<DataType> vector<DataType>::operator*(const double rhs) {
  return matrix<DataType>::operator*(rhs);
}

template<typename DataType>
vector<DataType>& vector<DataType>::operator*=(const double rhs) {
  return matrix<DataType>::operator*=(rhs);
}

template<typename DataType>
const vector<DataType> vector<DataType>::operator*(const CD rhs) {
  return matrix<DataType>::operator*(rhs);
}

template<typename DataType>
vector<DataType>& vector<DataType>::operator*=(const CD rhs) {
  return matrix<DataType>::operator*=(rhs);
}

template<typename DataType>
const vector<DataType> 
vector<DataType>::operator*(const matrix<DataType>& rhs) const {
  matrix<DataType> m = (*this);
  m.trans();
  vector<DataType> v1 = m*rhs;
  return v1;
}

template<typename DataType>
vector<DataType>& vector<DataType>::operator*=(const matrix<DataType>& rhs) {
  matrix<DataType> m = (*this);
  m.trans();
  (*this) = m*rhs;
  return *this;
}

template<typename DataType>
template<typename T> 
const DataType vector<DataType>::operator* (const vector<T> &rv) {
  assert(size() == rv.size());
  DataType res = 0;
  for(size_t i =0;i<size(); i++) res+=(*this)[i]*rv[i];
  return res;
}

template<typename DataType>
template<typename T> 
const vector<DataType> vector<DataType>::operator/(const T rhs) const {
  return matrix<DataType>::operator/(rhs);
}

template<typename DataType>
template<typename T> 
vector<DataType>& vector<DataType>::operator/=(const T rhs) {
  return matrix<DataType>::operator/=(rhs);
}

template<typename DataType>
vector<DataType> & vector<DataType>::trans() {return *this;};

template<typename DataType>
vector<DataType> & vector<DataType>::herm() {return conj();};

template<typename DataType>
vector<DataType> & vector<DataType>::conj() {
  matrix<DataType>::conj();
  return (*this);
}

template<typename DataType>
vector<DataType> vector<DataType>::block(size_t i1, size_t i2) {
  matrix<DataType>::block(i1, i2, 0, 1);
  return (*this);
}

template<typename DataType>
vector<DataType> & vector<DataType>::insert(size_t i, vector<DataType> v) {
  matrix<DataType>::insert(i, 0);
  return (*this);
}

template<typename DataType>
vector<DataType> & vector<DataType>::sort() {
  std::vector<DataType> v(matrix<DataType>::getDataPtr(), 
                          matrix<DataType>::getDataPtr()+size());
  std::sort(v.begin(), v.end());
  for (size_t i=0; i<v.size(); i++) (*this)[i] = v[i];
  return (*this);
}

template<typename DataType>
std::vector<DataType> vector<DataType>::stl() {
  std::vector<DataType> v(matrix<DataType>::getDataPtr(), 
                          matrix<DataType>::getDataPtr()+size());
  return v;
}

/////////////////////////////////////////////////////////////////////////////

template<typename DataType> 
vector<DataType> trans(const vector<DataType> &v) {
  return v;
}

template<typename DataType> 
vector<DataType> herm(const vector<DataType> &v) {
  vector<DataType> v1 = v;
  v1.herm();
  return v1;
}

template<typename DataType> 
vector<DataType> conj(const vector<DataType> &v) {
  vector<DataType> v1 = v;
  v1.conj();
  return v1;
}

vector<double> real(const vector<CD> &ma) {
 vector<double> m(ma.size());
 for (size_t i=0; i<ma.size(); i++) 
     m(i)=ma(i).real;
 return m;
} 

vector<double> imag(const vector<CD> &ma) {
 vector<double> m(ma.size());
 for (size_t i=0; i<ma.size(); i++) 
     m(i)=ma(i).imag;
 return m;
}

template<typename DataType> 
vector<DataType> sort(const vector<DataType> &v) {
  vector<DataType> v2(v.size());
  std::vector<DataType> v1(v.getDataPtr(), v.getDataPtr()+v.size());
  std::sort(v1.begin(), v1.end());
  for (size_t i=0; i<v1.size(); i++) v2[i] = v1[i];
  return v2;
}

template<typename DataType> 
matrix<DataType> diag(const vector<DataType> &v) {
  size_t n = v.size();
  matrix<DataType> m(n,n);
  for (size_t i=0; i<n; i++) 
    for (size_t j=0; j<n; j++)
      if (i==j) m(i,j)=v(i);
      else m(i,j)=0.;
  return m;
}

template<typename DataType> 
const vector<DataType> 
operator*(const double lhs, const vector<DataType> &ma) {
  vector<DataType> m(ma.size());
  for (size_t i=0; i<ma.nRow; i++)
      m(i) = ma(i)*lhs; 
  return m;
}

template<typename DataType> 
const vector<DataType> 
operator*(const CD lhs, const vector<DataType> &ma) {
  vector<DataType> m(ma.size());
  for (size_t i=0; i<ma.nRow; i++)
      m(i) = ma(i)*lhs; 
  return m;
}

template<typename DataType> 
const vector<DataType> 
operator*(const matrix<DataType>& ma, const vector<DataType>& v) {
  matrix<DataType> m = v;
  vector<DataType> v1 = ma*m;
  return v1;
}

template<typename DataType> 
vector<DataType>& 
operator*=(matrix<DataType>& ma, const vector<DataType>& v) {
  matrix<DataType> m = v;
  ma *= m;
  return ma;
}

template<typename DataType> 
std::ostream& operator<< (std::ostream& out, const vector<DataType> &m) {
  out << "{";
  for (size_t i=0; i<m.size(); i++) {
    if (i!=m.size()-1) out << m(i) << ", ";
    else out << m(i);
  }
  out << "}";
  return out; 
}

}
/////////////////////////////////////////////////////////////////////////////

#endif
