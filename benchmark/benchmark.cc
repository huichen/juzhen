/*
+---------------------------------------------------------------------------+
|  Juzhen: C++ library for linear algebra                                           |
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

int MAXN, MINN, STEPN;
char PREFIX[100]; 

int main(int argc, char *argv[]) {

  if (argc != 5) {
    MAXN = 100;
    MINN = 2;
    STEPN = 1;
    strcpy(PREFIX, "small_");
  } else {
    MAXN = atoi(argv[1]);
    MINN = atoi(argv[2]);
    STEPN = atoi(argv[3]);
    strcpy(PREFIX, argv[4]);
  }
  
  DEFV(1)
  DEFV(2)
  DEFV(3)
  DEFV(4)

  struct timeval t1, t2;
  FILE *myfile;
  char filename[100];

  BEGIN_BENCH("M1 = M2")
  M1 = M2;
  END_BENCH

  BEGIN_BENCH("M(i,j) = s, 10M times")
  for (size_t j = 0; j < 100000; j++) M1(ni/2, ni/2) = 100.;
  END_BENCH

  BEGIN_BENCH("M1 = M2.block (upper-left quarter), 1K times")
#ifdef MLCPP
  for (size_t j = 0; j < 100; j++) M1 = M2.Block(0, 0, ni/2-1, ni/2-1);
#else
  for (size_t j = 0; j < 100; j++) M1 = M2.block(0, 0, ni/2, ni/2);
#endif
  END_BENCH

  BEGIN_BENCH("M1 = M2.col(i), 1K times")
#ifdef MLCPP
  for (size_t j = 0; j < 100; j++) M1 = M2.GetCol(ni/2);
#else
  for (size_t j = 0; j < 100; j++) M1 = M2.col(ni/2);
#endif
  END_BENCH

  BEGIN_BENCH("M1 = M2.row(i), 1K times")
#ifdef MLCPP
  for (size_t j = 0; j < 100; j++) M1 = M2.GetRow(ni/2);
#else
  for (size_t j = 0; j < 100; j++) M1 = M2.row(ni/2);
#endif
  END_BENCH

  BEGIN_BENCH("M1 = M2.Conjugate()")
#ifdef MLCPP
  M1 = M2.Conjugate();
#else
  M1 = M2.conjugate();
#endif
  END_BENCH

  BEGIN_BENCH("M1 = M2.Transpose()")
#ifdef MLCPP
  M1 = M2.Transpose();
#else
  M1 = M2.transpose();
#endif
  END_BENCH

  BEGIN_BENCH("M1 = M2.Adjoint()")
#ifdef MLCPP
  M1 = M2.Adjoint();
#else
  M1 = M2.adjoint();
#endif
  END_BENCH

  BEGIN_BENCH("M1 = s*M2")
  M1 = 3.0*M2;
  END_BENCH

  BEGIN_BENCH("M1 = M2*s")
  M1 = M2*3.0;
  END_BENCH

  BEGIN_BENCH("M1 = M2; M1 *= s")
  M1 = M2;
  M1 *= 3.0;
  END_BENCH

  BEGIN_BENCH("M1 = M2/s")
  M1 = M2/3.0;
  END_BENCH

  BEGIN_BENCH("M1 = M2; M1 /= s")
  M1 = M2;
  M1 /= 3.0;
  END_BENCH

  BEGIN_BENCH("M1 = M2 + M3")
  M1 = M2 + M3;
  END_BENCH

  BEGIN_BENCH("M1 = s1*M2 + s2*M3")
  M1 = 2*M2 + 3*M3;
  END_BENCH

  BEGIN_BENCH("M1 = M2*M3")
  M1 = M2*M3;
  END_BENCH

  BEGIN_BENCH("M1 = s1*M2*M3 + s2*M4")
  M1 = 3*M2*M3 + 2*M4;
  END_BENCH

  return 0;
}
