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

#ifndef SRC_ADAPTOR_NATIVE_H_
#define SRC_ADAPTOR_NATIVE_H_
#include <omp.h>

#include <assert.h>

#include <adaptor/native/gemm.h>

namespace juzhen {

template<typename T>
int geev(
    char nl, char nr, const int n,
    T *a, const int lda, CD * w, T *vl, const int ldvl,
    T *vr, const int ldvr) {
  assert(0);  // always fails
}

/*
 * Linear solver
 */
template<typename T>
int gesv(
    const int n, const int nrhs,
    T *a, const int lda, T *b, const int ldb) {
  assert(0);  // always fails
}

/*
 * Matrix inversion 
 */
template<typename T>
int matrix_inverse(
    const int m, const int n, T *a, const int lda) {
  assert(0);  // always fails
}

/*
 * Matrix determinant 
 */
template<typename T>
T matrix_determinant(const int m, T *a) {
  assert(0);  // always fails
}
}
#endif  // SRC_ADAPTOR_NATIVE_H_
