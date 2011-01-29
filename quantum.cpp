#include "mklpp.h"
#define N 9
using namespace mklpp;

int main() {
  dvector potential(N);
  for (int i=0; i<N; i++) 
    if (i<N/3 || i>=2*N/3) potential(i) = 1000;
    else potential(i) = 0;

  cmatrix H = cdiagmatrix(potential);
  for (int i=0; i<N-1; i++) H(i, i+1) = H(i+1, i) = -1;

  cvector energy;
  cmatrix wave;
  H.reigen(energy, wave);

  std::cout << "energy= " <<  energy << std::endl;
  return 0;
}
