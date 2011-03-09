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
#include <juzhen.h>

#define N 90
using juzhen::cvector;
using juzhen::cmatrix;

int main() {
  /* matrix of complex<float> */
  cmatrix H(N, N);
  for (int i = 0; i < N; i++)
    if (i < N/3 || i >= 2*N/3)
      H(i, i) = 1000;
    else
      H(i, i) = 0;
  for (int i = 0; i < N-1; i++) H(i, i+1) = H(i+1, i) = -1;

#ifdef USE_MKL
  /* vector of complex<float> */
  cvector energy;
  /* matrix of complex<float> */
  cmatrix wave;
  RightEigenSolver(H, energy, wave);

  std::cout << "energy= " << Sort(Real(energy)) << std::endl;
#endif
  return 0;
}
