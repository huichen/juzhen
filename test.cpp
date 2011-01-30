#include "mklpp.h"
#include <vector>

#define P(s,m) std::cout << s << "=" << std::endl << m << std::endl;

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

  MKLCD c1[]= { {1.,0.}, {-3., 0.}, {2.0, 0.}, {-2., .0}, {6., 0.}, {3.0, 0.}};
  cmatrix m7(c1, 2,3);
  P("m7", m7)

  MKLCD c2[]={{ 1., 1.0},  { 2., 0.},    {3., 0.0 },   {6., 0.}, 
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

  MKLCD c3[] = {{1, 0}, {0, 0}, {0, 3.2}, {0, 0},
	{1.23, 0}, {1, 0}, {0, 0}, {0, 2},
	{3, 0.77}, {0, 0}, {1, 2.0}, {0, 0},
	{0, 4}, {2.2, 0}, {0, 0}, {1, 0}
	};
  cmatrix m10 (c3, 4, 4);
  m10.trans();
  P("m10", m10)
  cvector e;
  cmatrix vr ;
  m10.reigen(e, vr);
  P("eigen(m10):e", e)
  P("eigen(m10):vr", vr)
  P("diag(e)", diag(e))
  cmatrix ee = diag(e);

  P("ee", ee)
  P("m10*vr", m10*vr)
  P("vr*ee", vr*ee)
//  P("m10*vr-vr*e", m10*vr - vr*ee)

/* vector */
  double cv1[] = {1, 2, 3, 4};
  double cv2[] = {1, 3, 2, 4};
  
  dvector v1 (cv1, 4);
  dvector v2 (cv2, 4);

  P("v1", v1)
  P("v2", v2)
  P("v1*v1", v1*v1)
  P("v1+v2", v1+v2)
  P("3*v1", 3*v1)
  P("v1*3", v1*3)
  P("diag(v2)", diag(v2))
  P("diag(v2)*v1", diag(v2)*v1)
  P("v1*diag(v2)", v1*diag(v2))

  vector<double> sv1(4);
  sv1[0]=1;
  sv1[1]=2;
  sv1[2]=3;
  sv1[3]=7;

  P("dvector(sv1)", dvector(sv1))
  dvector v3;
  v3 = sv1;
  P("v3", v3)
  P("dvector(v3.stl())", dvector(v3.stl()))

  return 0;
}
