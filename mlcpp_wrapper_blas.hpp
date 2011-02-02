#ifndef MLCPP_WRAPPER_BLAS_HPP
#define MLCPP_WRAPPER_BLAS_HPP
#include <cblas.h>
#include <assert.h>

namespace mlcpp {

/* MKL function wrappers implemented with templates*/
template<typename T> 
void gemm(
  const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA, 
  const CBLAS_TRANSPOSE TransB, const int M, const int N, 
  const int K, const T *A, const int lda, const T *B, 
  const int ldb, T *c, const int ldc) { 
  assert(0); // always fails
};

template<>
void gemm<double>(
  const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA, 
  const CBLAS_TRANSPOSE TransB, const int M, const int N, 
  const int K, const double *A, const int lda, const double *B, 
  const int ldb, double *c, const int ldc) { 
  cblas_dgemm(Order, TransA, TransB, M, N, K, 1., 
              A, lda, B, ldb, 0., c, ldc); 
};

template<> 
void gemm<CD>(
  const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA, 
  const CBLAS_TRANSPOSE TransB, const int M, const int N, 
  const int K, const CD  *A, const int lda, const CD  *B, 
  const int ldb, CD  *c, const int ldc) { 
  CD alpha (1., 0.);
  CD beta (0., 0.);
  cblas_zgemm(Order, TransA, TransB, M, N, K, &alpha, 
              A, lda, B, ldb, &beta, c, ldc); 
};

template<typename T> 
int geev(
  const int order, char nl, char nr, const int n, 
  T *a, const int lda, CD * w, T *vl, const int ldvl, 
  T *vr, const int ldvr) { 
  assert(0); // always fails
};

}

#endif
