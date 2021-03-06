/*
+---------------------------------------------------------------------------+
|  Juzhen: C++ library for linear algebra                                   |
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

#ifndef SRC_CORE_COMPLEX_H_
#define SRC_CORE_COMPLEX_H_
#include <math.h>

#include <iostream>  // NOLINT string fuctions are provided for compatibility

namespace juzhen {

/////////////////////////////////////////////////////////////////////////////
/**
 * complex supplements the native complex structure with some
 * missing interfaces.
 */
template<typename T>
struct Complex {
  /**
   * Default constructor.
   */
  Complex() {}

  /**
   * Construct a Complex number from its real and imaginary part.
   */
  Complex(const T &r, const T &i);

  /**
   * Construct a Complex number from a real number with imaginary part to be
   * zero.
   */
  Complex(const T &r);  // NOLINT

  /**
   * Construct(copying) a Complex number from another Complex number.
   */
  template<typename T1>
  Complex(const Complex<T1> &c);  // NOLINT

  /**
   * Destructor.
   */
  ~Complex() {}

  /**
   * Copy from a real number with imaginary part set to zero.
   */
  inline const Complex<T> &operator=(const T &r);

  /**
   * Copy from another Complex number
   */
  template<typename T1>
  inline const Complex<T> &operator=(const Complex<T1> &c);

  /**
   * Return the complex number itself 
   */
  inline Complex<T> &operator+();

  /**
   * Return the opposite of a complex number 
   */
  inline Complex<T> operator-();

  /**
   * Check if two Complex numbers equal.
   */
  template<typename T1>
  inline bool operator==(const Complex<T1> &c);

  /**
   * Check if a Complex number equals to a real number and its imaginary part
   * is zero.
   */
  inline bool operator==(const T &c);

  /**
   * Check if two Complex numbers are not equal.
   */
  template<typename T1>
  inline bool operator!=(const Complex<T1> &c);

  /**
   * Check if a Complex number is not equal to a real number. It holds true
   * when the imaginary part of the Complex number is not zero.
   */
  inline bool operator!=(const T &c);

  /**
   * Real part of the complex number.
   */
  T real;

  /**
   * Imaginary part of the complex number.
   */
  T imag;
};

typedef Complex<double> CD;
typedef Complex<float> CS;

template<typename T>
Complex<T>::Complex(const T &r, const T &i) {
  real = r;
  imag = i;
}

template<typename T>
Complex<T>::Complex(const T &r) {
  real = r;
  imag = 0.0;
}

template<typename T>
template<typename T1>
Complex<T>::Complex(const Complex<T1> &c) {
  real = c.real;
  imag = c.imag;
}

template<typename T>
inline const Complex<T>& Complex<T>::operator=(const T &r) {
  real = r;
  imag = 0.0;
  return *this;
}

template<typename T>
template<typename T1>
inline const Complex<T> &Complex<T>::operator=(const Complex<T1> &c) {
  real = c.real;
  imag = c.imag;
  return *this;
}

template<typename T>
inline Complex<T> &Complex<T>::operator+() {
  return *this;
}

template<typename T>
inline Complex<T> Complex<T>::operator-() {
  return Complex<T>(-real, -imag);
}

template<typename T>
template<typename T1>
inline bool Complex<T>::operator==(const Complex<T1> &c) {
  if (real == c.real && imag == c.imag)
    return true;
  else
    return false;
}

template<typename T>
inline bool Complex<T>::operator==(const T &c) {
  if (real == c && imag == 0)
    return true;
  else
    return false;
}

template<typename T>
template<typename T1>
inline bool Complex<T>::operator!=(const Complex<T1> &c) {
  if (real == c.real && imag == c.imag)
    return false;
  else
    return true;
}

template<typename T>
inline bool Complex<T>::operator!=(const T &c) {
  if (real == c && imag == 0)
    return false;
  else
    return true;
}



/////////////////////////////////////////////////////////////////////////////

/* arithmatic operations */
template<typename T>
inline const Complex<T> operator+(const Complex<T> &a, const Complex<T> &b) {
  Complex<T> c;
  c.real = a.real + b.real;
  c.imag = a.imag + b.imag;
  return c;
}

template<typename T>
inline const Complex<T> operator+(const Complex<T> &a, const T &b) {
  Complex<T> c;
  c.real = a.real + b;
  return c;
}

template<typename T>
inline Complex<T> &operator+=(Complex<T> &a, const Complex<T> &b) {
  a.real += b.real;
  a.imag += b.imag;
  return a;
}

template<typename T>
inline Complex<T>& operator+=(Complex<T> &a, const T &b) {
  a.real += b;
  return a;
}

template<typename T>
inline const Complex<T> operator-(const Complex<T> &a, const Complex<T> &b) {
  Complex<T> c;
  c.real = a.real - b.real;
  c.imag = a.imag - b.imag;
  return c;
}

template<typename T>
inline const Complex<T> operator-(const Complex<T> &a, const T &b) {
  Complex<T> c;
  c.real = a.real - b;
  return c;
}

template<typename T>
inline Complex<T>& operator-=(Complex<T> &a, const Complex<T> &b) {
  a.real -= b.real;
  a.imag -= b.real;
  return a;
}

template<typename T>
inline Complex<T>& operator-=(Complex<T> &a, const T &b) {
  a.real -= b;
  return a;
}

template<typename T>
inline const Complex<T> operator*(const Complex<T> &a, const Complex<T> &b) {
  Complex<T> c;
  c.real = a.real * b.real - a.imag * b.imag;
  c.imag = a.real * b.imag + a.imag * b.real;
  return c;
}

template<typename T, typename T1>
inline Complex<T> &operator*=(Complex<T> &a, const Complex<T1> &b) {
  a.real = a.real * b.real - a.imag * b.imag;
  a.imag = a.real * b.imag + a.imag * b.real;
  return a;
}

template<typename T, typename T1>
inline const Complex<T> operator*(const Complex<T> &a, const T1 &b) {
  Complex<T> c;
  c.real = a.real * b;
  c.imag = a.imag * b;
  return c;
}

template<typename T, typename T1>
inline const Complex<T> operator*(const T1 &b, const Complex<T> &a) {
  Complex<T> c;
  c.real = a.real * b;
  c.imag = a.imag * b;
  return c;
}

template<typename T, typename T1>
inline Complex<T> &operator*=(Complex<T> &a, const T1 &b) {
  a.real = a.real * b;
  a.imag = a.imag * b;
  return a;
}

template<typename T, typename T1>
inline const Complex<T> operator/(const Complex<T> &a, const T1 &b) {
  Complex<T> c;
  c.real = a.real / b;
  c.imag = a.imag / b;
  return c;
}

template<typename T, typename T1>
inline Complex<T>& operator/=(Complex<T> &a, const T1 &b) {
  a.real = a.real / b;
  a.imag = a.imag / b;
  return a;
}

template<typename T, typename T1>
inline const Complex<T> operator/(const T1 &b, const Complex<T> &a) {
  Complex<T> c;
  T temp = b / (a.real * a.real + a.imag * a.imag);
  c.real = a.real * temp;
  c.imag = -a.imag * temp;
  return c;
}

template<typename T>
inline T Abs2(const Complex<T> &a) {
  return a.real * a.real + a.imag * a.imag;
}

inline double Abs2(double a) {
  return a * a;
}

inline float Abs2(float a) {
  return a * a;
}

template<typename T>
inline T Abs(const Complex<T> &a) {
  return sqrt(a.real * a.real + a.imag * a.imag);
}

inline double Abs(double a) {
  return a > 0 ? a : -a;
}

inline float Abs(float a) {
  return a > 0 ? a : -a;
}

template<typename T, typename T1>
inline const Complex<T> operator/(const Complex<T> &a, const Complex<T1> &b) {
  Complex<T> c;
  T de = 1 / Abs2(b);
  c.real = de * (a.real * b.real + a.imag * b.imag);
  c.imag = de * (-a.real * b.imag + a.imag * b.real);
  return c;
}

template<typename T, typename T1>
inline Complex<T> &operator/=(Complex<T> &a, const Complex<T1> &b) {
  T de = 1 / Abs2(b);
  T temp = a.real;
  a.real = de * (a.real * b.real + a.imag * b.imag);
  a.imag = de * (-temp * b.imag + a.imag * b.real);
  return a;
}

template<typename T>
inline const Complex<T> Conjugate(const Complex<T> &a) {
  Complex<T> c;
  c.real = a.real;
  c.imag = -a.imag;
  return c;
}

inline double Conjugate(double a) {
  return a;
}

inline float Conjugate(float a) {
  return a;
}

template<typename T>
inline T Real(const Complex<T> &a) {
  return a.real;
}

template<typename T>
inline T Imag(const Complex<T> &a) {
  return a.imag;
}

inline float Real(float a) {
  return a;
}

template<typename T>
inline float Imag(float a) {
  return a;
}

inline double Real(double a) {
  return a;
}

template<typename T>
inline double Imag(double a) {
  return a;
}

template<typename T>
std::ostream &operator<< (std::ostream &out, const Complex<T> &m) {
  out << "(" << m.real << ", " << m.imag << ")";
  return out;
}
/////////////////////////////////////////////////////////////////////////////
}
#endif  // SRC_CORE_COMPLEX_H_

