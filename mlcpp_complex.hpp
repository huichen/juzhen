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

#ifdef USE_MKL
#include <mkl_cblas.h>
#include <mkl_lapacke.h>
#include <mkl.h>
#else
#include <cblas.h>
#endif

namespace mlcpp {

#ifdef USE_MKL
#define _CD MKL_Complex16 
#else
/**
 * _CD is the native complex struct that contains two double numbers, real
 * part before the imaginary part.
 */
typedef struct { 
  /**
   * Real part of the complex number.
   */
  double real; 
  /**
   * Imaginary part of the complex number.
   */
  double imag;
} _CD; 
#endif

#define CD ComplexDouble 

/////////////////////////////////////////////////////////////////////////////
/** 
 * _ComplexDouble supplements the native complex structure with some
 * missing interfaces.
 */
struct _ComplexDouble : _CD {

  /** 
   * Default constructor.
   */
  _ComplexDouble() {};

  /** 
   * Construct a complex number from its real and imaginary part. 
   */
  _ComplexDouble(double r, double i);

  /** 
   * Construct a complex number from a real number with imaginary part to be
   * zero.
   */
  _ComplexDouble(double r);

  /**
   * Construct(copying) a complex number from another complex number.
   */
  _ComplexDouble(_CD c);

  /**
   * Copy from a real number with imaginary part set to zero. 
   */
  inline const struct _ComplexDouble & operator=(double r);

  /**
   * Copy from another complex number
   */
  inline const struct _ComplexDouble & operator=(const _CD &c);

  /**
   * Check if two complex numbers equal.
   */
  inline bool operator==(const _CD &c);

  /**
   * Check if a complex number equals to a real number and its imaginary part
   * is zero. 
   */
  inline bool operator==(const double c);

  /**
   * Check if two complex numbers are not equal.
   */
  inline bool operator!=(const _CD &c);

  /**
   * Check if a complex number is not equal to a real number. It holds true
   * when the imaginary part of the complex number is not zero. 
   */
  inline bool operator!=(const double c);

};

typedef struct _ComplexDouble CD;

_ComplexDouble::_ComplexDouble(double r, double i) { real = r; imag = i;};
_ComplexDouble::_ComplexDouble(double r) { real = r; imag =0.0;};
_ComplexDouble::_ComplexDouble(_CD c) { real = c.real; imag =c.imag;};

inline const CD& _ComplexDouble::operator=(double r) {
  real = r;
  imag = 0.0;
  return *this;
}

inline const CD& _ComplexDouble::operator=(const _CD &c) {
  real = c.real;
  imag = c.imag;
  return *this;
}

inline bool _ComplexDouble::operator==(const _CD &c) {
  if (real == c.real && imag == c.imag) return true;
  else return false;
}

inline bool _ComplexDouble::operator==(const double c) {
  if (real == c && imag == 0) return true;
  else return false;
}

inline bool _ComplexDouble::operator!=(const _CD &c) {
  if (real == c.real && imag == c.imag) return false;
  else return true;
}

inline bool _ComplexDouble::operator!=(const double c) {
  if (real == c && imag == 0) return false;
  else return true;
}



/////////////////////////////////////////////////////////////////////////////

/* complex constant */
CD C1 (1.,0.);
CD Ci (0.,1.);
CD C0 (0.,0.);

/* arithmatic operations */
inline const CD operator+(const CD &a, const CD &b) {
  CD c;
  c.real = a.real + b.real; 
  c.imag = a.imag + b.imag;
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

inline const CD operator*(const double b, const CD &a) {
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

inline const CD operator/(const CD &a, const double b) {
  CD c;
  c.real = a.real / b; 
  c.imag = a.imag / b;
  return c;
}

inline CD& operator/=(CD &a, const double b) {
  a.real = a.real / b; 
  a.imag = a.imag / b;
  return a;
}

inline double abs2(const CD &a) {
  return a.real*a.real+a.imag*a.imag;
}

inline const CD operator/(const CD &a, const CD &b) {
  CD c;
  double de = 1/abs2(b);
  c.real = de*(a.real*b.real + a.imag*b.imag);
  c.imag = de*(-a.real*b.imag + a.imag*b.real);
  return c;
}

inline CD& operator/=(CD &a, const CD &b) {
  double de = 1/abs2(b);
  double temp = a.real;
  a.real = de*(a.real*b.real + a.imag*b.imag);
  a.imag = de*(-temp*b.imag + a.imag*b.real);
  return a;
}

inline const CD conj(const CD &a) {
  CD c;
  c.real = a.real; 
  c.imag = -a.imag;
  return c;
}

inline double conj(const double &a) {
  return a;
}

std::ostream& operator<< (std::ostream& out, const CD &m) {
  out << "(" << m.real << ", " << m.imag << ")";
  return out; 
}
/////////////////////////////////////////////////////////////////////////////
 
}
#endif
