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

#ifndef MLCPP_ADAPTOR_BLAS_H_  // NOLINT
#define MLCPP_ADAPTOR_BLAS_H_
#include <assert.h>

#include <cblas.h>

namespace mlcpp {

/* MKL function wrappers implemented with templates*/
template<typename T>
void gemm(
    const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
    const CBLAS_TRANSPOSE TransB, const int M, const int N,
    const int K, const T *A, const int lda, const T *B,
    const int ldb, T *c, const int ldc) {
  assert(0);  // always fails
}

template<>
void gemm<float>(
    const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
    const CBLAS_TRANSPOSE TransB, const int M, const int N,
    const int K, const float *A, const int lda, const float *B,
    const int ldb, float *c, const int ldc) {
  cblas_sgemm(Order, TransA, TransB, M, N, K, 1.,
              A, lda, B, ldb, 0., c, ldc);
}

template<>
void gemm<double>(
    const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
    const CBLAS_TRANSPOSE TransB, const int M, const int N,
    const int K, const double *A, const int lda, const double *B,
    const int ldb, double *c, const int ldc) {
  cblas_dgemm(Order, TransA, TransB, M, N, K, 1.,
              A, lda, B, ldb, 0., c, ldc);
}

template<>
void gemm<CS>(
    const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
    const CBLAS_TRANSPOSE TransB, const int M, const int N,
    const int K, const CS  *A, const int lda, const CS *B,
    const int ldb, CS  *c, const int ldc) {
  CS alpha(1., 0.);
  CS beta(0., 0.);
  cblas_cgemm(Order, TransA, TransB, M, N, K, &alpha,
              A, lda, B, ldb, &beta, c, ldc);
}

template<>
void gemm<CD>(
    const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
    const CBLAS_TRANSPOSE TransB, const int M, const int N,
    const int K, const CD  *A, const int lda, const CD  *B,
    const int ldb, CD  *c, const int ldc) {
  CD alpha(1., 0.);
  CD beta(0., 0.);
  cblas_zgemm(Order, TransA, TransB, M, N, K, &alpha,
              A, lda, B, ldb, &beta, c, ldc);
}

template<typename T>
int geev(
    const int order, char nl, char nr, const int n,
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
#endif  // MLCPP_ADAPTOR_BLAS_H_  // NOLINT

