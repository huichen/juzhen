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

#ifndef MLCPP_ADAPTOR_MKL_H_  // NOLINT
#define MLCPP_ADAPTOR_MKL_H_
#include <assert.h>

#include <mkl_cblas.h>
#include <mkl_lapacke.h>
#include <mkl.h>

namespace mlcpp {

/* MKL function wrappers implemented with templates*/
template<typename T>
void gemm(
    const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
    const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N,
    const MKL_INT K, const T *A, const MKL_INT lda, const T *B,
    const MKL_INT ldb, T *c, const MKL_INT ldc) {
  assert(0);  // always fails
}

template<>
void gemm<float>(
    const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
    const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N,
    const MKL_INT K, const float *A, const MKL_INT lda, const float *B,
    const MKL_INT ldb, float *c, const MKL_INT ldc) {
  cblas_sgemm(Order, TransA, TransB, M, N, K, 1.,
              A, lda, B, ldb, 0., c, ldc);
}

template<>
void gemm<double>(
    const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
    const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N,
    const MKL_INT K, const double *A, const MKL_INT lda, const double *B,
    const MKL_INT ldb, double *c, const MKL_INT ldc) {
  cblas_dgemm(Order, TransA, TransB, M, N, K, 1.,
              A, lda, B, ldb, 0., c, ldc);
}

template<>
void gemm<CS>(
    const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
    const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N,
    const MKL_INT K, const CS  *A, const MKL_INT lda, const CS  *B,
    const MKL_INT ldb, CS  *c, const MKL_INT ldc) {
  CS alpha(1., 0.);
  CS beta(0., 0.);
  cblas_cgemm(Order, TransA, TransB, M, N, K, &alpha,
              A, lda, B, ldb, &beta, c, ldc);
}

template<>
void gemm<CD>(
    const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
    const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N,
    const MKL_INT K, const CD  *A, const MKL_INT lda, const CD  *B,
    const MKL_INT ldb, CD  *c, const MKL_INT ldc) {
  CD alpha(1., 0.);
  CD beta(0., 0.);
  cblas_zgemm(Order, TransA, TransB, M, N, K, &alpha,
              A, lda, B, ldb, &beta, c, ldc);
}

template<typename T, typename T1>
int geev(
    char nl, char nr, const MKL_INT n,
    T *a, const MKL_INT lda, T1 * w, T *vl, const MKL_INT ldvl,
    T *vr, const MKL_INT ldvr) {
  assert(0);  // always fails
}

template<>
int geev<float, CS>(
    char nl, char nr, const MKL_INT n,
    float *a, const MKL_INT lda, CS * w, float *vl, const MKL_INT ldvl,
    float *vr, const MKL_INT ldvr) {
  float *wr = reinterpret_cast<float *>(malloc(sizeof(float)*n));  // NOLINT
  float *wi = reinterpret_cast<float *>(malloc(sizeof(float)*n));  // NOLINT
  assert(wr);
  assert(wi);
  int res = LAPACKE_sgeev(LAPACK_COL_MAJOR, nl, nr, n,
                         a, lda, wr, wi, vl, ldvl, vr, ldvr);
  for (size_t i = 0; i < n; i++) {w[i].real = wr[i]; w[i].imag = wi[i];}
  free(wr);
  free(wi);
  return res;
}

template<>
int geev<double, CD>(
    char nl, char nr, const MKL_INT n,
    double *a, const MKL_INT lda, CD * w, double *vl, const MKL_INT ldvl,
    double *vr, const MKL_INT ldvr) {
  double *wr = reinterpret_cast<double *>(malloc(sizeof(double)*n));  // NOLINT
  double *wi = reinterpret_cast<double *>(malloc(sizeof(double)*n));  // NOLINT
  assert(wr);
  assert(wi);
  int res = LAPACKE_dgeev(LAPACK_COL_MAJOR, nl, nr, n,
                         a, lda, wr, wi, vl, ldvl, vr, ldvr);
  for (size_t i = 0; i < n; i++) {w[i].real = wr[i]; w[i].imag = wi[i];}
  free(wr);
  free(wi);
  return res;
}

template<>
int geev<CS, CS>(
    char nl, char nr, const MKL_INT n,
    CS *a, const MKL_INT lda, CS * w, CS *vl, const MKL_INT ldvl,
    CS *vr, const MKL_INT ldvr) {
  return LAPACKE_cgeev(
      LAPACK_COL_MAJOR, nl, nr, n,
      reinterpret_cast<MKL_Complex8 *>(a), lda,
      reinterpret_cast<MKL_Complex8 *>(w),
      reinterpret_cast<MKL_Complex8 *>(vl), ldvl,
      reinterpret_cast<MKL_Complex8 *>(vr), ldvr);
}

template<>
int geev<CD, CD>(
    char nl, char nr, const MKL_INT n,
    CD *a, const MKL_INT lda, CD * w, CD *vl, const MKL_INT ldvl,
    CD *vr, const MKL_INT ldvr) {
  return LAPACKE_zgeev(
      LAPACK_COL_MAJOR, nl, nr, n,
      reinterpret_cast<MKL_Complex16 *>(a), lda,
      reinterpret_cast<MKL_Complex16 *>(w),
      reinterpret_cast<MKL_Complex16 *>(vl), ldvl,
      reinterpret_cast<MKL_Complex16 *>(vr), ldvr);
}
}

#endif  // MLCPP_ADAPTOR_MKL_H_  // NOLINT
