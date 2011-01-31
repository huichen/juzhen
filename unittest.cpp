#include <iostream>
#include <mklpp.h>
#include <unittest.hpp>
#include <vector>
#include <sstream>

using namespace mklpp;

INIT_UNITTEST

int main() {

BEGIN_TEST(MatrixAssignment, "MatrixAssignment")
  double a[] = {1., 2., 3., 4.};
  double b[] = {1., 2., 3., 4., 5., 6.};

  matrix<double> m1(a, 2, 2);
  V(m1, "1 3 \n2 4 ")

  matrix<double> m2 = m1;
  V(m2, "1 3 \n2 4 ")

  m2(0,1) = 99.;
  V(m2, "1 99 \n2 4 ")

  m1 = m2;
  V(m1, "1 99 \n2 4 ")

  matrix<double> m3(b, 3, 2);
  V(m3, "1 4 \n2 5 \n3 6 ")

  matrix<double> m4(3,3);
  m4 = 3;
  V(m4, "3 3 3 \n3 3 3 \n3 3 3 ")

  cmatrix m5 = m3;
  V(m5, "(1, 0) (4, 0) \n(2, 0) (5, 0) \n(3, 0) (6, 0) ")

  cmatrix m6;
  m6 = m3;
  V(m6, "(1, 0) (4, 0) \n(2, 0) (5, 0) \n(3, 0) (6, 0) ")

  cmatrix m7(m6.getDataPtr(),2,2);
  V(m7, "(1, 0) (3, 0) \n(2, 0) (4, 0) ")

END_TEST(MatrixAssignment)

BEGIN_TEST(IdentityMatrix, "IdentityMatrix")
  matrix<double> m1 = idmatrix<double>(3);
  V(m1, "1 0 0 \n0 1 0 \n0 0 1 ")

  cidmatrix m2(2);
  V(m2, "(1, 0) (0, 0) \n(0, 0) (1, 0) ")

END_TEST(IdentityMatrix)

BEGIN_TEST(MatrixHermitian, "MatrixHermitian")

  MKLCD c1[]= { {1.,0.}, {-3., 0.}, {2.0, 0.}, {-2., .0}, {6., 0.}, {3.0, 0.}};
  cmatrix m1(c1, 2,3);
  V(herm(m1), "(1, 0) (-3, 0) \n(2, 0) (-2, 0) \n(6, 0) (3, 0) ")

  VE(m1, herm(herm(m1)))

END_TEST(MatrixHermitian)

/*
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
*/
BEGIN_TEST(VectorOperation, "VectorOperation")
  double cv1[] = {1, 2, 3, 4};
  double cv2[] = {1, 3, 2, 4};
  
  dvector v1 (cv1, 4);
  dvector v2 (cv2, 4);
  
/*
  P("v1*v1", v1*v1)
  P("v1+v2", v1+v2)
  P("3*v1", 3*v1)
  P("v1*3", v1*3)
  P("diag(v2)", diag(v2))
  P("diag(v2)*v1", diag(v2)*v1)
  P("v1*diag(v2)", v1*diag(v2))
*/
  vector<double> sv1(4);
  sv1[0]=1;
  sv1[1]=2;
  sv1[2]=3;
  sv1[3]=7;

  V(sv1,"{1, 2, 3, 7}")
  V(dvector(sv1),"{1, 2, 3, 7}")

  dvector v3;
  v3 = sv1;
  V(v3,"{1, 2, 3, 7}")

  V(dvector(v3.stl()), "{1, 2, 3, 7}")

END_TEST(VectorOperation)

  UnitTest MyTest;
  MyTest.Run(); 

  return 0;
}