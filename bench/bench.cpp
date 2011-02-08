#include "bench.hpp"

int main() {

  DEFV(1)
  DEFV(2)
  DEFV(3)
  DEFV(4)

  struct timeval t1, t2;
  std::ofstream myfile;
  std::string filename;

  BEGIN_BENCH ("M1 = M2")
    H1 = H2; 
  END_BENCH

  BEGIN_BENCH ("M(i,j)")
  for (size_t j=0; j<1000000; j++) H1(3,4) = H2(5,6) + H3(7,8);
  END_BENCH

  BEGIN_BENCH ("M1 = M2.block")
  for (size_t j=0; j<100; j++) H1 = H2.block(0,0,50,50);
  END_BENCH

  BEGIN_BENCH ("M1 = M2.col(i)")
  for (size_t j=0; j<1000; j++)H1 = H2.col(10);
  END_BENCH

  BEGIN_BENCH ("M1 = M2.row(i)")
  H1 = H2.row(10);
  END_BENCH

  BEGIN_BENCH ("M1 = M2.conj()")
  H1 = conj(H2);
  END_BENCH

  BEGIN_BENCH ("M1 = M2.trans()")
  #ifdef MLCPP
  H1 = trans(H2);
  #else
  H1 = H2.transpose();
  #endif
  END_BENCH

  BEGIN_BENCH ("M1 = M2.herm()")
  #ifdef MLCPP
  H1 = herm(H2);
  #else
  H1 = H2.adjoint();
  #endif
  END_BENCH

  BEGIN_BENCH ("M1 = s*M2")
  H1 = 3.0*H2;
  END_BENCH

  BEGIN_BENCH ("M1 = M2*s")
  H1 = H2*3.0;
  END_BENCH

  BEGIN_BENCH ("M1 = M2; M1 *= s")
  H1 = H2;
  H1 *= 3.0;
  END_BENCH

  BEGIN_BENCH ("M1 = M2/s")
  H1 = H2/3.0;
  END_BENCH

  BEGIN_BENCH ("M1 = M2; M1 /= s") 
  H1 = H2;
  H1 /= 3.0;
  END_BENCH

  BEGIN_BENCH ("M1 = M2 + M3")
  H1 = H2 + H3;
  END_BENCH

  BEGIN_BENCH ("M1 = s1*M2 + s2*M3")
  H1 = 2*H2 + 3*H3;
  END_BENCH

/*
  BEGIN_BENCH ("M1 = M2*M3")
  H1 = H2*H3;
  END_BENCH

  BEGIN_BENCH ("M1 = s1*M2*M3 + s2*M4")
  H1 = 3*H2*H3 + 2*H4;
  END_BENCH
*/
  return 0;
}
