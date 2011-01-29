#include "mklpp.h"

#define P(s,m) cout << s << "=" << endl << m << endl;

using namespace mklpp;

int main() {

// assign and copy

  double a[] = {1., 2., 3., 4.};
  double b[] = {1., 2., 3., 4., 5., 6.};

  matrix<double> m1(a, 2, 2);
  P("m1(0,1)", m1(0,1))

  matrix<double> m2 = m1;
  m2(0,1) = 99.;

  matrix<double> m3(b, 2, 3);

  P("m1",m1)
  P("m2",m2)

  matrix<double> m4(3,3);
  m4 = 3;
  P("m4",m4)

// identity matrix
  matrix<double> m5 = idmatrix<double>(3);
  P("idmatrix(3)=", idmatrix<double>(3))
  P("m5", m5)
  P("m4*m5", m4*m5)

  cmatrix m6(4,4);
  m6 = C1*3. + Ci*2.;

  P("m6", m6)
  P("1", cidmatrix(4))
  P("m6*1", m6*cidmatrix(4))

// complex conjugate transpose

  CD c1[]={ {1.,0.}, {-3., 0.}, {2.0, 0.}, {-2., .0}, {6., 0.}, {3.0, 0.}};
  cmatrix m7(c1, 2,3);
  P("m7", m7)

  CD c2[]={{ 1., 1.0},  { 2., 0.},    {3., 0.0 },   {6., 0.}, 
           { 1., .0},  { 2., 0.},    {3., 0.0 },   {6., 0.}, 
           { 1., .0},  { 2., 0.},    {3., 0.0 },   {6., 0.}};
  cmatrix m8(c2, 4, 3);
  m8.herm();
  P("m8", m8)

// arithmetic operations

  P("m7*m8", m7*m8)
  P("m7*m8", m7*m8)
  P("m7+m7", m7+m7)
  
  cmatrix m9 = m7;
  m9+=m7;
  P("m9=m7; (m9+=m7)", m9)
  m9-=m7;
  P("(m9-=m7)", m9)

  P("2*m7.", 2.*m7)
  P("m7/4.", m7/4.)

// sub-matrix
  P("m7(1->2, 1->3)", m7.block(1,2,1,3))
  P("m8(1->2, 1->2)=I ", m8.insert(1,1,cidmatrix(2)))
  P("m7(1,)", m7.row(1))
  P("m7(,1)", m7.col(1))

// eigen systems

  CD c3[] = {{1, 0}, {0, 0}, {0, 3.2}, {0, 0},
	{1.23, 0}, {1, 0}, {0, 0}, {0, 2},
	{3, 0.77}, {0, 0}, {1, 2.0}, {0, 0},
	{0, 4}, {2.2, 0}, {0, 0}, {1, 0}
	};
  cmatrix m10 (c3, 4, 4);
  m10.trans();
  P("m10", m10)
  cmatrix e, vr ;
  m10.reigen(e, vr);
  P("eigen(m10):e", e)
  P("eigen(m10):vr", vr)

  cmatrix ee = cdiagmatrix(e);

  P("ee", ee)
  P("m10*vr", m10*vr)
  P("vr*ee", vr*ee)
//  P("m10*vr-vr*e", m10*vr - vr*ee)
  return 0;
}
