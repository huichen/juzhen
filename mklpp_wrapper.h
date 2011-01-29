#ifndef MKLPP_WRAPPER_H
#define MKLPP_WRAPPER_H
#include <mkl_cblas.h>
#include <mkl_lapacke.h>
#include <mkl.h>

namespace mklpp {

/* MKL function wrappers implemented with templates*/
template<typename T> void gemm(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N, const MKL_INT K, const T *A, const MKL_INT lda, const T *B, const MKL_INT ldb, T *c, const MKL_INT ldc) { 
  assert(0); // always fails
};

template<> void gemm<double>(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N, const MKL_INT K, const double *A, const MKL_INT lda, const double *B, const MKL_INT ldb, double *c, const MKL_INT ldc) { 
  cblas_dgemm(Order, TransA, TransB, M, N, K, 1., A, lda, B, ldb, 0., c, ldc); 
};

template<> void gemm<CD>(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N, const MKL_INT K, const CD  *A, const MKL_INT lda, const CD  *B, const MKL_INT ldb, CD  *c, const MKL_INT ldc) { 
  CD alpha (1., 0.);
  CD beta (0., 0.);
  cblas_zgemm(Order, TransA, TransB, M, N, K, &alpha, A, lda, B, ldb, &beta, c, ldc); 
};

template<typename T> int geev(const int order, char nl, char nr, const MKL_INT n, T *a, const MKL_INT lda, T * w, T *vl, const MKL_INT ldvl, T *vr, const MKL_INT ldvr) { 
  assert(0); // always fails
};

template<> int geev<CD>(const int order, char nl, char nr, const MKL_INT n, CD *a, const MKL_INT lda, CD * w, CD *vl, const MKL_INT ldvl, CD *vr, const MKL_INT ldvr) { 
  return LAPACKE_zgeev(LAPACK_COL_MAJOR, nl, nr, n, a, lda, w, vl, ldvl, vr, ldvr); 
};


}

#endif
