#include <mlcpp.hpp>
#define N 90
using namespace mlcpp;

int main() {
  dvector potential(N);
  for (int i=0; i<N; i++) 
    if (i<N/3 || i>=2*N/3) potential(i) = 1000;
    else potential(i) = 0;

  cmatrix H = diag(potential);
  for (int i=0; i<N-1; i++) H(i, i+1) = H(i+1, i) = -1;

#ifdef USE_MKL
  cvector energy;
  cmatrix wave;
  H.reigen(energy, wave);

  std::cout << "energy= " << sort(real(energy)) << std::endl;
#endif
  return 0;
}
