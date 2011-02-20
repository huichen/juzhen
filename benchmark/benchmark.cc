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
#include "benchmark.h"  // NOLINT

int main() {
  DEFV(1)
  DEFV(2)
  DEFV(3)
  DEFV(4)

  struct timeval t1, t2;
  FILE *myfile;
  char filename[100];

  BEGIN_BENCH("M1 = M2")
    H1 = H2;
  END_BENCH

  BEGIN_BENCH("M(i,j), 1M times")
  for (size_t j = 0; j < 100000; j++) H1(3, 4) = H2(5, 6) + H3(7, 8);
  END_BENCH

  BEGIN_BENCH("M1 = M2.block, 1K times")
#ifdef MLCPP
  for (size_t j = 0; j < 100; j++) H1 = H2.Block(0, 0, 50, 50);
#else
  for (size_t j = 0; j < 100; j++) H1 = H2.block(0, 0, 50, 50);
#endif
  END_BENCH

  BEGIN_BENCH("M1 = M2.col(i), 1K times")
#ifdef MLCPP
  for (size_t j = 0; j < 100; j++) H1 = H2.GetCol(1);
#else
  for (size_t j = 0; j < 100; j++) H1 = H2.col(1);
#endif
  END_BENCH

  BEGIN_BENCH("M1 = M2.row(i), 1K times")
#ifdef MLCPP
  for (size_t j = 0; j < 100; j++) H1 = H2.GetRow(1);
#else
  for (size_t j = 0; j < 100; j++) H1 = H2.row(1);
#endif
  END_BENCH

  BEGIN_BENCH("M1 = M2.Conjugate()")
#ifdef MLCPP
  H1 = Conjugate(H2);
#else
  H1 = H2.conjugate();
#endif
  END_BENCH

  BEGIN_BENCH("M1 = M2.Transpose()")
#ifdef MLCPP
  H1 = Transpose(H2);
#else
  H1 = H2.transpose();
#endif
  END_BENCH

  BEGIN_BENCH("M1 = M2.Adjoint()")
#ifdef MLCPP
  H1 = Adjoint(H2);
#else
  H1 = H2.adjoint();
#endif
  END_BENCH

  BEGIN_BENCH("M1 = s*M2")
  H1 = 3.0*H2;
  END_BENCH

  BEGIN_BENCH("M1 = M2*s")
  H1 = H2*3.0;
  END_BENCH

  BEGIN_BENCH("M1 = M2; M1 *= s")
  H1 = H2;
  H1 *= 3.0;
  END_BENCH

  BEGIN_BENCH("M1 = M2/s")
  H1 = H2/3.0;
  END_BENCH

  BEGIN_BENCH("M1 = M2; M1 /= s")
  H1 = H2;
  H1 /= 3.0;
  END_BENCH

  BEGIN_BENCH("M1 = M2 + M3")
  H1 = H2 + H3;
  END_BENCH

  BEGIN_BENCH("M1 = s1*M2 + s2*M3")
  H1 = 2*H2 + 3*H3;
  END_BENCH

  BEGIN_BENCH("M1 = M2*M3")
  H1 = H2*H3;
  END_BENCH

  BEGIN_BENCH("M1 = s1*M2*M3 + s2*M4")
  H1 = 3*H2*H3 + 2*H4;
  END_BENCH

  return 0;
}
