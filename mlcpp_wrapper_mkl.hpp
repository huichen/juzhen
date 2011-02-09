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

#ifndef MLCPP_WRAPPER_MKL_HPP
#define MLCPP_WRAPPER_MKL_HPP
#include <mkl_cblas.h>
#include <mkl_lapacke.h>
#include <mkl.h>
#include <assert.h>

namespace mlcpp {

/* MKL function wrappers implemented with templates*/
template<typename T> 
void gemm(
  const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA, 
  const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N, 
  const MKL_INT K, const T *A, const MKL_INT lda, const T *B, 
  const MKL_INT ldb, T *c, const MKL_INT ldc) { 
  assert(0); // always fails
};

template<>
void gemm<double>(
  const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA, 
  const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N, 
  const MKL_INT K, const double *A, const MKL_INT lda, const double *B, 
  const MKL_INT ldb, double *c, const MKL_INT ldc) { 
  cblas_dgemm(Order, TransA, TransB, M, N, K, 1., 
              A, lda, B, ldb, 0., c, ldc); 
};

template<> 
void gemm<CD>(
  const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA, 
  const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N, 
  const MKL_INT K, const CD  *A, const MKL_INT lda, const CD  *B, 
  const MKL_INT ldb, CD  *c, const MKL_INT ldc) { 
  CD alpha (1., 0.);
  CD beta (0., 0.);
  cblas_zgemm(Order, TransA, TransB, M, N, K, &alpha, 
              A, lda, B, ldb, &beta, c, ldc); 
};

template<typename T> 
int geev(
  char nl, char nr, const MKL_INT n, 
  T *a, const MKL_INT lda, CD * w, T *vl, const MKL_INT ldvl, 
  T *vr, const MKL_INT ldvr) { 
  assert(0); // always fails
};

template<> 
int geev<CD>(
  char nl, char nr, const MKL_INT n, 
  CD *a, const MKL_INT lda, CD * w, CD *vl, const MKL_INT ldvl, 
  CD *vr, const MKL_INT ldvr) { 
  return LAPACKE_zgeev(LAPACK_COL_MAJOR, nl, nr, n, 
                       a, lda, w, vl, ldvl, vr, ldvr); 
};

template<> 
int geev<double>(
  char nl, char nr, const MKL_INT n, 
  double *a, const MKL_INT lda, CD * w, double *vl, const MKL_INT ldvl, 
  double *vr, const MKL_INT ldvr) { 

  double *wr = (double*) malloc(sizeof(double)*n); 
  double *wi = (double*) malloc(sizeof(double)*n); 
  assert(wr);
  assert(wi);
  int res = LAPACKE_dgeev(LAPACK_COL_MAJOR, nl, nr, n, 
                         a, lda, wr, wi, vl, ldvl, vr, ldvr); 
  for(size_t i=0; i<n; i++) {w[i].real = wr[i]; w[i].imag = wi[i];}
  free(wr);
  free(wi);
  return res;
};

}

#endif
