#ifndef MLCPP_VECTOR_HPP
#define MLCPP_VECTOR_HPP

#include <mlcpp_matrix.hpp>
#include <algorithm>
#include <vector>
#include <math.h>

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
  inline size_t size() const;

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
  const DataType operator* (const vector<T> &rv);

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
   * Get real part of a vector.
   */
  vector<double> real() const;

  /**
   * Get imaginary part of a vector.
   */
  vector<double> imag() const;

  /**
   * Get transpose of a vector (just conjugate it). 
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
   * Get sum of a vector's all elements.
   */
  DataType sum() const; 

  /**
   * Get norm_2 of a vector.
   */
  double norm() const; 

  /**
   * Get the maximum element of a vector.
   */
  DataType max() const; 

  /**
   * Swap two elements.
   */
  vector<DataType> &swap(size_t i1, size_t i2);

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
                                                    m.nrow()*m.ncol())); 
  matrix<DataType>::m_ncol =1;
  matrix<DataType>::m_nrow = m.nrow()*m.ncol();
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
  matrix<DataType>::operator=(v);
  return *this;
} 

template<typename DataType>
vector<DataType> & 
vector<DataType>::operator=(const std::vector<DataType> &v) { 
  size_t s = v.size();
  resize(s);
  DataType *p1;
  p1 = matrix<DataType>::getDataPtr();
  for(size_t i=0; i<s; i++) *(p1++) = v[i];
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
inline size_t vector<DataType>::size() const {
  return matrix<DataType>::m_nrow;
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
const DataType vector<DataType>::operator* (const vector<T> &rv) {
  assert(size() == rv.size());
  DataType res = 0;
  size_t s = size();
  DataType *p1 = matrix<DataType>::getDataPtr();
  T *p2 = rv.getDataPtr();
  for(size_t i =0;i<s; i++) res+=*(p1++)*(*(p2++));
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
vector<double> vector<DataType>::real() const {
  return matrix<DataType>::real();
}

template<typename DataType>
vector<double> vector<DataType>::imag() const {
  return matrix<DataType>::imag();
}

template<typename DataType>
vector<DataType> vector<DataType>::trans() const {
  return *this;
}


template<typename DataType>
DataType vector<DataType>::sum() const {
  DataType r = 0;
  DataType *p = matrix<DataType>::getDataPtr();
  for(size_t i=0; i<size(); i++)
    r += *(p++);
  return r;
} 

template<typename DataType>
double vector<DataType>::norm() const {
  double r = 0;
  DataType *p = matrix<DataType>::getDataPtr();
  for(size_t i=0; i<size(); i++)
    r += abs2(*(p++));
  return sqrt(r);
} 

template<typename DataType>
DataType vector<DataType>::max() const {
  assert(size()>0);
  DataType r = (*this)(0);
  DataType *p = matrix<DataType>::getDataPtr();
  for(size_t i=0; i<size(); i++) {
    if (r < *p) r = *p; 
    p++;
  }
  return r;
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
vector<DataType> & vector<DataType>::swap(size_t i1, size_t i2) {
  assert(i1<size() && i2<size());
  DataType *p = matrix<DataType>::m_data->m_data;
  DataType temp;
  temp = p[i1];
  p[i1] = p[i2];
  p[i2] = temp;
  return *this;
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

/**
 * Return the transpose of a vector.
 */
template<typename DataType> 
vector<DataType> trans(const vector<DataType> &v) {
  return v;
}

/**
 * Return the hermitian of a vector.
 */
template<typename DataType> 
vector<DataType> herm(const vector<DataType> &v) {
  return v.herm();
}

/**
 * Return the conjugate of a vector.
 */
template<typename DataType> 
vector<DataType> conj(const vector<DataType> &v) {
  return v.conj();
}

/**
 * Return the real part of a vector.
 */
vector<double> real(const vector<CD> &v) {
  return v.real();
} 

/**
 * Return the imaginary part of a vector.
 */
vector<double> imag(const vector<CD> &v) {
  return v.imag();
}

/**
 * Find the sum of a vector's all elements.
 */
template<typename DataType> 
DataType sum(const vector<DataType> v) {
  return v.sum();
} 

/**
 * Find the norm-2 of a vector.
 */
template<typename DataType> 
double norm(const vector<DataType> v) {
  return v.norm();
} 

/**
 * Find the maximum element in a vector.
 */
template<typename DataType> 
DataType max(const vector<DataType> v) {
  return v.max();
} 

/**
 * Sort a vector.
 */
template<typename DataType> 
vector<DataType> sort(const vector<DataType> &v) {
  vector<DataType> v2(v.size());
  std::vector<DataType> v1(v.getDataPtr(), v.getDataPtr()+v.size());
  std::sort(v1.begin(), v1.end());
  for (size_t i=0; i<v1.size(); i++) v2[i] = v1[i];
  return v2;
}

/**
 * Multiply a real number and a vector.
 */
template<typename DataType> 
const vector<DataType> 
operator*(const double lhs, const vector<DataType> &ma) {
  vector<DataType> m(ma.size());
  for (size_t i=0; i<ma.size(); i++)
    m(i) = ma(i)*lhs; 
  return m;
}

/**
 * Multiply a complex number and a vector.
 */
template<typename DataType> 
const vector<DataType> 
operator*(const CD lhs, const vector<DataType> &ma) {
  vector<DataType> m(ma.size());
  for (size_t i=0; i<ma.size(); i++)
      m(i) = ma(i)*lhs; 
  return m;
}

/**
 * Multiply a matrix and a vector.
 */
template<typename DataType> 
const vector<DataType> 
operator*(const matrix<DataType>& ma, const vector<DataType>& v) {
  matrix<DataType> m = v;
  vector<DataType> v1 = ma*m;
  return v1;
}

/**
 * Multiply a matrix and a vector.
 */
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

/**
 * Return string form of a matrix.
 */

template<typename DataType> 
std::string toString(const vector<DataType> &m) {
  std::ostringstream out;
  out << m;
  return out.str(); 
}


}
/////////////////////////////////////////////////////////////////////////////

#endif
