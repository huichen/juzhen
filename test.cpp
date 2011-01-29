#include "mkl++.h"

#define P(s,m) cout << s << "=" << endl << m << endl;

using namespace MAT;

int main() {

// assign and copy

  double a[] = {1., 2., 3., 4.};
  double b[] = {1., 2., 3., 4., 5., 6.};

  Matrix<double> m1(a, 2, 2);
  P("m1(0,1)", m1(0,1))

  Matrix<double> m2 = m1;
  m2(0,1) = 99.;

  Matrix<double> m3(b, 2, 3);

  P("m1",m1)
  P("m2",m2)

  Matrix<double> m4(3,3);
  m4 = 3;
  P("m4",m4)

// identity matrix
  Matrix<double> m5 = IdentityMatrix<double>(3);
  cout << m5.nCol << m5.nRow << sizeof(m5.getDataPtr()) << endl;
  P("IdentityMatrix(3)=", IdentityMatrix<double>(3))
  P("m5", m5)
  P("m4*m5", m4*m5)

  CMatrix m6(4,4);
  m6 = C1*3. + Ci*2.;

  P("m6", m6)
  P("1", CIdentityMatrix(4))
  P("m6*1", m6*CIdentityMatrix(4))

// complex conjugate transpose

  CD c1[]={ CD(1.,0.), CD(-3., 0.), CD(2.0, 0.), CD(-2., .0), CD(6., 0.), CD(3.0, 0.)};
  CMatrix m7(c1, 2,3);
  P("m7", m7)

  CD c2[]={CD( 1., 1.0),  CD( 2., 0.),    CD(3., 0.0 ),   CD(6., 0.), 
           CD( 1., .0),  CD( 2., 0.),    CD(3., 0.0 ),   CD(6., 0.), 
           CD( 1., .0),  CD( 2., 0.),    CD(3., 0.0 ),   CD(6., 0.)};
  CMatrix m8(c2, 4, 3);
  m8.conjTrans();
  P("m8", m8)

// arithmetic operations

  P("m7*m8", m7*m8)
  P("m7*m8", m7*m8)
  P("m7+m7", m7+m7)
  
  CMatrix m9 = m7;
  m9+=m7;
  P("m9=m7; (m9+=m7)", m9)
  m9-=m7;
  P("(m9-=m7)", m9)

  P("2*m7.", 2.*m7)
  P("m7/4.", m7/4.)

// sub-matrix
  P("m7(1->2, 1->3)", m7.subMatrix(1,2,1,3))
  P("m8(1->2, 1->2)=I ", m8.replace(1,1,CIdentityMatrix(2)))
  P("m7(1,)", m7.rowVector(1))
  P("m7(,1)", m7.colVector(1))

  return 0;
}
