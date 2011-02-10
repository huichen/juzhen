/*
+---------------------------------------------------------------------------+
|  Matrix Library for C++ (mlcpp)                                           |
+---------------------------------------------------------------------------+
|                                                                           |
|  Copyright 2011 Hui Chen                                                  |
|                                                                           |
|  Licensed under the Apache License, Version 2.0 (the "License");          |
|  you may not use this file except in compliance with the License.         |
|  You may obtain a copy of the License at                                  |
|                                                                           |
|      http://www.apache.org/licenses/LICENSE-2.0                           |
|                                                                           |
|  Unless required by applicable law or agreed to in writing, software      |
|  distributed under the License is distributed on an "AS IS" BASIS,        |
|  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. |
|  See the License for the specific language governing permissions and      |
|  limitations under the License.                                           |
|                                                                           |
+---------------------------------------------------------------------------+
*/

#ifndef MLCPP_COMPLEX_HPP
#define MLCPP_COMPLEX_HPP
#include <iostream>

namespace mlcpp {

/////////////////////////////////////////////////////////////////////////////
/** 
 * complex supplements the native complex structure with some
 * missing interfaces.
 */
template<typename T>
struct complex {
  /**
   * Real part of the complex number.
   */
  T real; 
  /**
   * Imaginary part of the complex number.
   */
  T imag;

  /** 
   * Default constructor.
   */
  complex() {};

  /** 
   * Construct a complex number from its real and imaginary part. 
   */
  complex(T r, T i);

  /** 
   * Construct a complex number from a real number with imaginary part to be
   * zero.
   */
  complex(const T r);

  /**
   * Construct(copying) a complex number from another complex number.
   */
  template<typename T1>
  complex(const complex<T1>& c);

  /**
   * Copy from a real number with imaginary part set to zero. 
   */
  inline const struct complex & operator=(T r);

  /**
   * Copy from another complex number
   */
  template<typename T1>
  inline const struct complex & operator=(const complex<T1> &c);

  /**
   * Check if two complex numbers equal.
   */
  template<typename T1>
  inline bool operator==(const complex<T1> &c);

  /**
   * Check if a complex number equals to a real number and its imaginary part
   * is zero. 
   */
  inline bool operator==(const T c);

  /**
   * Check if two complex numbers are not equal.
   */
  template<typename T1>
  inline bool operator!=(const complex<T1> &c);

  /**
   * Check if a complex number is not equal to a real number. It holds true
   * when the imaginary part of the complex number is not zero. 
   */
  inline bool operator!=(const T c);

  typedef T ValueType;
};

typedef complex<double> CD;
typedef complex<float> CS;

template<typename T>
complex<T>::complex(T r, T i) { 
  real = r; 
  imag = i;
}

template<typename T>
complex<T>::complex(T r) { 
  real = r; 
  imag =0.0;
}

template<typename T>
template<typename T1>
complex<T>::complex(const complex<T1>& c) { 
  real = c.real; 
  imag =c.imag;
}

template<typename T>
inline const complex<T>& complex<T>::operator=(T r) {
  real = r;
  imag = 0.0;
  return *this;
}

template<typename T>
template<typename T1>
inline const complex<T>& complex<T>::operator=(const complex<T1> &c) {
  real = c.real;
  imag = c.imag;
  return *this;
}

template<typename T>
template<typename T1>
inline bool complex<T>::operator==(const complex<T1> &c) {
  if (real == c.real && imag == c.imag) return true;
  else return false;
}

template<typename T>
inline bool complex<T>::operator==(const T c) {
  if (real == c && imag == 0) return true;
  else return false;
}

template<typename T>
template<typename T1>
inline bool complex<T>::operator!=(const complex<T1> &c) {
  if (real == c.real && imag == c.imag) return false;
  else return true;
}

template<typename T>
inline bool complex<T>::operator!=(const T c) {
  if (real == c && imag == 0) return false;
  else return true;
}



/////////////////////////////////////////////////////////////////////////////

/* arithmatic operations */
template<typename T>
inline const complex<T> operator+(const complex<T> &a, const complex<T> &b) {
  complex<T> c;
  c.real = a.real + b.real; 
  c.imag = a.imag + b.imag;
  return c;
}

template<typename T>
inline complex<T>& operator+=(complex<T> &a, const complex<T> &b) {
  a.real += b.real; 
  a.imag += b.imag;
  return a;
}

template<typename T>
inline const complex<T> operator-(const complex<T> &a, const complex<T> &b) {
  complex<T> c;
  c.real = a.real - b.real; 
  c.imag = a.imag - b.imag;
  return c;
}

template<typename T>
inline complex<T>& operator-=(complex<T> &a, const complex<T> &b) {
  a.real -= b.real; 
  a.imag -= b.real;
  return a;
}

template<typename T>
inline const complex<T> operator*(const complex<T> &a, const complex<T> &b) {
  complex<T> c;
  c.real = a.real*b.real - a.imag*b.imag; 
  c.imag = a.real*b.imag + a.imag*b.real;
  return c;
}

template<typename T, typename T1>
inline complex<T>& operator*=(complex<T> &a, const complex<T1> &b) {
  a.real = a.real*b.real - a.imag*b.imag; 
  a.imag = a.real*b.imag + a.imag*b.real;
  return a;
}

template<typename T, typename T1>
inline const complex<T> operator*(const complex<T> &a, const T1 b) {
  complex<T> c;
  c.real = a.real * b; 
  c.imag = a.imag * b;
  return c;
}

template<typename T, typename T1>
inline const complex<T> operator*(const T1 b, const complex<T> &a) {
  complex<T> c;
  c.real = a.real * b; 
  c.imag = a.imag * b;
  return c;
}

template<typename T, typename T1>
inline complex<T>& operator*=(complex<T> &a, const T1 b) {
  a.real = a.real * b; 
  a.imag = a.imag * b;
  return a;
}

template<typename T, typename T1>
inline const complex<T> operator/(const complex<T> &a, const T1 b) {
  complex<T> c;
  c.real = a.real / b; 
  c.imag = a.imag / b;
  return c;
}

template<typename T, typename T1>
inline complex<T>& operator/=(complex<T> &a, const T1 b) {
  a.real = a.real / b; 
  a.imag = a.imag / b;
  return a;
}

template<typename T>
inline T abs2(const complex<T> &a) {
  return a.real*a.real+a.imag*a.imag;
}

inline double abs2(const double a) {
  return a*a;
}

inline float abs2(const float a) {
  return a*a;
}

template<typename T, typename T1>
inline const complex<T> operator/(const complex<T> &a, const complex<T1> &b) {
  complex<T> c;
  T de = 1/abs2(b);
  c.real = de*(a.real*b.real + a.imag*b.imag);
  c.imag = de*(-a.real*b.imag + a.imag*b.real);
  return c;
}

template<typename T, typename T1>
inline complex<T>& operator/=(complex<T> &a, const complex<T1> &b) {
  T de = 1/abs2(b);
  T temp = a.real;
  a.real = de*(a.real*b.real + a.imag*b.imag);
  a.imag = de*(-temp*b.imag + a.imag*b.real);
  return a;
}

template<typename T>
inline const complex<T> conj(const complex<T> &a) {
  complex<T> c;
  c.real = a.real; 
  c.imag = -a.imag;
  return c;
}

inline double conj(const double &a) {
  return a;
}

inline float conj(const float &a) {
  return a;
}

template<typename T>
std::ostream& operator<< (std::ostream& out, const complex<T> &m) {
  out << "(" << m.real << ", " << m.imag << ")";
  return out; 
}
/////////////////////////////////////////////////////////////////////////////
 
}
#endif
