#ifndef MLCPP_VECTOR_HPP
#define MLCPP_VECTOR_HPP

#include <mlcpp_matrix.hpp>
#include <algorithm>
#include <vector>

namespace mlcpp {

/////////////////////////////////////////////////////////////////////////////
/** 
 *  vector class
 */
template<typename DataType> class vector : public matrix<DataType> {
public:
  /**
   * Default constructor.
   */
  vector();

  /**
   * Construct a vector of size i.
   */
  vector(size_t i);

  /**
   * Construct a vector from a matrix's column-major raw array.
   */
  vector(const matrix<DataType> &m);

  /**
   * Construct a vector from an array. The size of the vector will be
   * the bigger of data's size and n.
   */
  template<typename T> 
  vector(const T *data, size_t n);

  /**
   * Construct a vector from an existing vector.
   */
  vector(const std::vector<DataType> &v);

  /**
   * Copy a vector from another vector.
   */
  vector<DataType>& operator=(const vector<DataType> &v);

  /**
   * Copy a vector from a STL vector.
   */
  vector<DataType>& operator=(const std::vector<DataType> &v);

  /**
   * Assign all elements in the vector to be rhs. 
   */
  vector<DataType>& operator=(const DataType rhs);

  /**
   * Resize a vector.
   */
  void resize(size_t n);

  /**
   * This has no effect. 
   */
  void resize(size_t nc, size_t nr);

  /**
   * Return the size of the vector.
   */
  size_t size() const;

  /**
   * Add two vectors.
   */
  template<typename T> 
  const vector<DataType> operator+(const vector<T>& rhs) const;
 
  /**
   * Add two vectors.
   */
  template<typename T> 
  vector<DataType>& operator+=(const vector<T>& rhs);
 
  /**
   * Subtract two vectors.
   */
  template<typename T> 
  const vector<DataType> operator-(const vector<T>& rhs) const;
 
  /**
   * Subtract two vectors.
   */
  template<typename T> 
  vector<DataType>& operator-=(const vector<T>& rhs); 

  /**
   * Multiply a vector with a real number. 
   */
  const vector<DataType> operator*(const double rhs) const;

  /**
   * Multiply a vector with a real number. 
   */
  vector<DataType>& operator*=(const double rhs); 

  /**
   * Multiply a vector with a complex number. 
   */
  const vector<DataType> operator*(const CD rhs) const;

  /**
   * Multiply a vector with a complex number. 
   */
  vector<DataType>& operator*=(const CD rhs);

  /**
   * Multiply a vector and a matrix. The vector's size must be equal
   * to the matrix's number of rows. 
   */
  const vector<DataType> operator*(const matrix<DataType>& rhs) const;

  /**
   * Multiply a vector and a matrix. The vector's size must be equal
   * to the matrix's number of rows. 
   */
  vector<DataType>& operator*=(const matrix<DataType>& rhs);

  /**
   * Dot product two vectors. 
   */
  template<typename T> 
  const DataType operator* (const vector<T> &rv) const;

  /**
   * Divide a vector by a constant.
   */
  template<typename T> 
  const vector<DataType> operator/(const T rhs) const;

  /**
   * Divide a vector by a constant.
   */
  template<typename T> 
  vector<DataType>& operator/=(const T rhs); 

  /**
   * Get transpose of a vector (no effect).
   */
  vector<DataType> trans() const;

  /**
   * Get hermitian of a vector (just conjugate it). 
   */
  vector<DataType> herm() const;

  /**
   * Get conjugate of a vector.
   */
  vector<DataType> conj() const; 

  /**
   * Return a sub-vector of [i1, i2)
   */
  vector<DataType> block(size_t i1, size_t i2);

  /**
   * Replace the sub-vector starting at i with another vector v. 
   */
  vector<DataType> &insert(size_t i, vector<DataType> v);

  /**
   * Sort a vector. This will give you error for a complex vector.
   */
  vector<DataType> &sort();

  /**
   * Get a STL vector from the vector. 
   */
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
  matrix<DataType>::m_data = 
   typename DataPtr<DataType>::Type(new DataArray<DataType>(m.getDataPtr(), 
                                                    m.nRow()*m.nCol())); 
  matrix<DataType>::m_ncol =1;
  matrix<DataType>::m_nrow = m.nRow()*m.nCol();
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
  matrix<DataType>::m_data = 
   typename DataPtr<DataType>::Type(new DataArray<DataType>(v.getDataPtr(), 
                                                    v.size())); 
  matrix<DataType>::m_ncol =1;
  matrix<DataType>::m_nrow = v.size();
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
vector<DataType>& vector<DataType>::operator=(const DataType rhs) {
  matrix<DataType>::operator=(rhs);
  return *this;
}

template<typename DataType>
void vector<DataType>::resize(size_t n) { matrix<DataType>::resize(n,1);}

template<typename DataType>
void vector<DataType>::resize(size_t nc, size_t nr) {};

template<typename DataType>
size_t vector<DataType>::size() const {
  return matrix<DataType>::m_ncol==1?matrix<DataType>::m_nrow:matrix<DataType>::m_ncol;
}

template<typename DataType>
template<typename T> 
const vector<DataType> vector<DataType>::operator+(const vector<T>& rhs) const {
  return matrix<DataType>::operator+(rhs);
}

template<typename DataType>
template<typename T> 
vector<DataType>& vector<DataType>::operator+=(const vector<T>& rhs) {
  matrix<DataType>::operator+=(rhs);
  return *this;
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
  matrix<DataType>::operator-=(rhs);
  return *this;
}

template<typename DataType>
const vector<DataType> vector<DataType>::operator*(const double rhs) const {
  return matrix<DataType>::operator*(rhs);
}

template<typename DataType>
vector<DataType>& vector<DataType>::operator*=(const double rhs) {
  matrix<DataType>::operator*=(rhs);
  return *this;
}

template<typename DataType>
const vector<DataType> vector<DataType>::operator*(const CD rhs) const {
  return matrix<DataType>::operator*(rhs);
}

template<typename DataType>
vector<DataType>& vector<DataType>::operator*=(const CD rhs) {
  matrix<DataType>::operator*=(rhs);
  return *this;
}

template<typename DataType>
const vector<DataType> 
vector<DataType>::operator*(const matrix<DataType>& rhs) const {
  matrix<DataType> m = (*this);
  vector<DataType> v1 = m.trans() * rhs;
  return v1;
}

template<typename DataType>
vector<DataType>& vector<DataType>::operator*=(const matrix<DataType>& rhs) {
  matrix<DataType> m = (*this);
  (*this) = trans(m)*rhs;
  return *this;
}

template<typename DataType>
template<typename T> 
const DataType vector<DataType>::operator* (const vector<T> &rv) const {
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
  matrix<DataType>::operator/=(rhs);
  return *this;
}

template<typename DataType>
vector<DataType> vector<DataType>::trans() const {
  return *this;
}

template<typename DataType>
vector<DataType> vector<DataType>::herm() const {
  return matrix<DataType>::conj();
}

template<typename DataType>
vector<DataType> vector<DataType>::conj() const {
  return matrix<DataType>::conj();
}

template<typename DataType>
vector<DataType> vector<DataType>::block(size_t i1, size_t i2) {
  return matrix<DataType>::block(i1, i2, 0, 1);
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
  return v.herm();
}

template<typename DataType> 
vector<DataType> conj(const vector<DataType> &v) {
  return v.conj();
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
  for (size_t i=0; i<ma.nRow(); i++)
      m(i) = ma(i)*lhs; 
  return m;
}

template<typename DataType> 
const vector<DataType> 
operator*(const CD lhs, const vector<DataType> &ma) {
  vector<DataType> m(ma.size());
  for (size_t i=0; i<ma.nRow(); i++)
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

template<typename DataType> 
std::string toString(const vector<DataType> &m) {
  std::ostringstream out;
  out << m;
  return out.str(); 
}


}
/////////////////////////////////////////////////////////////////////////////

#endif
