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

#ifndef MLCPP_ADAPTOR_NATIVE_H_  // NOLINT
#define MLCPP_ADAPTOR_NATIVE_H_
#include <omp.h>

#include <assert.h>

namespace mlcpp {

template<typename T>
void gemm(
    const int M, const int N,
    const int K, const T *A, const int lda, const T *B,
    const int ldb, T *c, const int ldc) {
  const T *ai, *bi;
  T *ci;
  #pragma omp parallel for private(ai, bi, ci)
  for (size_t i = 0; i < M; i++) {
    for (size_t j = 0; j < N; j++) {
      ai = A + i;
      bi = B + j * ldb;
      ci = c + j * ldc + i;
      (*ci) = 0;
      for (size_t k = 0; k < K; k++) {
        (*ci) += (*ai) * (*(bi++));
        ai += lda;
      }
    }
  }
}

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
}
#endif  // MLCPP_ADAPTOR_NATIVE_H_  // NOLINT

