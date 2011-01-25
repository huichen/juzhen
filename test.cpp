#include "mkl++.h"

using namespace MAT;

int main() {

  Matrix<double> m1;

  double a[] = {1., 2., 3., 4.};
  double b[] = {1., 2., 3., 4., 5., 6.};

  m1 = a;
  m1.setDim(2, 2);
  std::cout << m1(0,1) << std::endl;

  Matrix<double> m2 = m1;
  m2(0,1) = 99.;

  Matrix<double> m3(b, 2, 3);

  cout << "m1=" << endl << m1 << endl;
  cout << "m2=" << endl << m2 << endl;
  cout << "m1*m2=" << endl <<m1*m2 << endl;

  cout << "m2=" << endl <<m3 << endl;
  cout << "m1*m1=" << endl <<m1*m1 << endl;

  Matrix<double> m4 = DMatrix<double>(3);

  cout << "m4=" << endl << m4 << endl;

  Matrix<double> m5(3,3);
  m5 = 3;

  cout << "m4*m5=" << endl << m4*m5 << endl;

  CMatrix m6(4,4);
  m6 = C1*3. + Ci*2.;

  cout << "m6=" << endl << m6 << endl;
  cout << "1=" << endl << CDMatrix(4) << endl;
  cout << "m6*1=" << endl << m6*CDMatrix(4) << endl;

  CD c1[]={ CD(1.,0.), CD(-3., 0.), CD(2.0, 0.), CD(-2., .0), CD(6., 0.), CD(3.0, 0.)};
  CMatrix m7(c1, 2,3);
  //cout << "t(m7)=" << endl << transpose(m7) << endl;
  //cout << "t^2(m7)=" << endl << transpose(transpose(m7)) << endl;
  cout << "m7=" << endl << m7 << endl;

  CD c2[]={CD( 1., 1.0),  CD( 2., 0.),    CD(3., 0.0 ),   CD(6., 0.), 
           CD( 1., .0),  CD( 2., 0.),    CD(3., 0.0 ),   CD(6., 0.), 
           CD( 1., .0),  CD( 2., 0.),    CD(3., 0.0 ),   CD(6., 0.)};
  CMatrix m8(c2, 4, 3);
  m8.conjTrans();
  cout << "m8=" << endl << m8 << endl;

  cout << "m7*m8=" << endl << m7*m8 << endl;
  return 0;
}
