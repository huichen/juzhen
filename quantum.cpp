#include "mklpp.h"

#define P(s,m) cout << s << "=" << endl << m << endl;
#define N 6

using namespace mklpp;

int main() {

  matrix<double> potential(N, 1);
  for (int i=0; i<N; i++) 
    if (i<N/3 || i>=2*N/3) potential(i) = 1000.0;
    else potential(i) = 0.0;

  cmatrix H = cdiagmatrix(potential);
  for (int i=0; i<N-1; i++) {
    H(i, i+1) = -1.;
    H(i+1, i) = -1.;
  }

  cmatrix energy, wave;
  H.reigen(energy, wave);

  P("H", H)
  P("energy", energy)
  P("wave", wave)
  return 0;
}
