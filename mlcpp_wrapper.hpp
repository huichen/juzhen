#ifndef MLCPP_WRAPPER
#define MLCPP_WRAPPER

#ifdef USE_MKL
#include <mlcpp_wrapper_mkl.hpp>
#else
#include <mlcpp_wrapper_blas.hpp>
#endif

#endif
