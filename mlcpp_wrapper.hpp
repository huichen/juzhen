#ifndef MLCPP_WRAPPER
#define MLCPP_WRAPPER

#ifdef USE_MKL
#include <mlcpp_wrapper_mkl.hpp>
#else
#ifdef USE_ATLAS
#include <mlcpp_wrapper_atlas.hpp>
#endif
#endif

#endif
