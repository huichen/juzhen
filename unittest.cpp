#include <iostream>
#include <mlcpp.hpp>
#include <unittest.hpp>
#include <vector>
#include <sstream>

using namespace mlcpp;

int main() {

/**
 * Tests for matrix
*/
BEGIN_TEST(MatrixAssignment, "MatrixAssignment")
  double a[] = {1., 2., 3., 4.};
  double b[] = {1., 2., 3., 4., 5., 6.};

  dmatrix m1(a, 2, 2);
  V(m1, "1 3 \n2 4 ")

  dmatrix m2 = m1;
  V(m2, "1 3 \n2 4 ")

  m2(0,1) = 99.;
  V(m2, "1 99 \n2 4 ")

  m1 = m2;
  V(m1, "1 99 \n2 4 ")

  dmatrix m3(b, 3, 2);
  V(m3, "1 4 \n2 5 \n3 6 ")

  dmatrix m4(3,3);
  m4 = 3;
  V(m4, "3 3 3 \n3 3 3 \n3 3 3 ")

  cmatrix m5 = m3;
  V(m5, "(1, 0) (4, 0) \n(2, 0) (5, 0) \n(3, 0) (6, 0) ")

  cmatrix m6;
  m6 = m3;
  V(m6, "(1, 0) (4, 0) \n(2, 0) (5, 0) \n(3, 0) (6, 0) ")

  cmatrix m7(m6.getDataPtr(),2,2);
  V(m7, "(1, 0) (3, 0) \n(2, 0) (4, 0) ")

  dmatrix m8(4,4);
  m8 = 3;
  V(m8, "3 3 3 3 \n3 3 3 3 \n3 3 3 3 \n3 3 3 3 ")

  m8 = m8;
  V(m8, "3 3 3 3 \n3 3 3 3 \n3 3 3 3 \n3 3 3 3 ")

END_TEST(MatrixAssignment)

BEGIN_TEST(MatrixCompare, "MatrixCompare")

  double a[] = {1., 2., 3., 4.};
  double b[] = {1., 2., 3., 4., 5., 6.};

  dmatrix m1(a, 2, 2);
  dmatrix m2(a, 2, 3);
  
  VB(m1==m1)

  VB(m1!=m2)

  dmatrix m3 = m1;
  VB(m1==m3)

  m3(1, 1) = 5;
  VB(m1!=m3)

  cmatrix m4 = m1;
/* 
  It's not comparable between cmatrix and dmatrix
  VB(m4!=m1);
*/

  cmatrix m5;
  m5 = m4;
  VB(m4==m5)

  m4(1,1) = (_CD){3,3};
  VB(m4!=m5)

END_TEST(MatrixCompare)

BEGIN_TEST(MatrixPrivateData, "MatrixPrivateData")
  double a[] = {1., 2., 3., 4., 5., 6.};
  dmatrix m1(a, 2, 3);

  VDE(m1.nCol(), 3)

  VDE(m1.nRow(), 2)

END_TEST(MatrixPrivateData)

BEGIN_TEST(MatrixResize, "MatrixResize")
  double a[] = {1., 2., 3., 4., 5., 6.};
  dmatrix m1(a, 2, 3);
 
  m1.resize(3,2); 
  V(m1, "1 4 \n2 5 \n3 6 ")
 
  m1.resize(3,1); 
  V(m1, "1 \n2 \n3 ")

  m1.resize(3,2); 
  V(m1, "1 4 \n2 5 \n3 6 ")

  m1.resize(10,10); 
  m1.resize(3,2); 
  V(m1, "1 4 \n2 5 \n3 6 ")

END_TEST(MatrixResize)

BEGIN_TEST(MatrixIndex, "MatrixIndex")
  double a[] = {1., 2., 3., 4., 5., 6.};
  dmatrix m1(a, 2, 3);

  VDE(m1(1, 1), 4) 

  m1(1, 2) = 7;
  VDE(m1(1, 2), 7) 

  VDE(m1(4), 5) 

  m1(2) = 9;
  VDE(m1(2), 9)

  VDE(m1(0, 1), 9)

END_TEST(MatrixIndex)

BEGIN_TEST(MatrixArithmatic, "MatrixArithmatic")
  double a[] = {1., 2., 3., 4., 5., 6.};
  dmatrix m1(a, 2, 3);

  dmatrix m2 = m1+m1;
  V(m2, "2 6 10 \n4 8 12 ")

  m2 = m1;
  m2 += m1;
  V(m2, "2 6 10 \n4 8 12 ")

  m2 = m1 + m1;
  m2 = m2 - m1; 
  V(m2, "1 3 5 \n2 4 6 ")

  m2 = m1 + m1;
  m2 -= m1; 
  V(m2, "1 3 5 \n2 4 6 ")

  V(m1*2, "2 6 10 \n4 8 12 ")

  V(2*m1, "2 6 10 \n4 8 12 ")

  m2 = m1;
  m2 *= 2;
  V(m2, "2 6 10 \n4 8 12 ")

  V(m2/2, "1 3 5 \n2 4 6 ")

  m2 /= 2;
  V(m2, "1 3 5 \n2 4 6 ")

  double b[] = {1, 2, 3, 4};
  dmatrix m3(b, 2, 2);

  dmatrix m4;
  m4 = m3*m2;
  V(m4, "7 15 23 \n10 22 34 ")

  m4 = m3;
  m4 *= m2;
  V(m4, "7 15 23 \n10 22 34 ")

END_TEST(MatrixArithmatic)

BEGIN_TEST(MatrixHermTransConj, "MatrixHermTransConj")

  _CD c1[]= { {1.,2.}, {-3., -1.}, {2.0, 0.1}, {-2., 0.1}, {6., 7.}, {3.0, 0.1}};
  cmatrix m1(c1, 2,3);
  V(m1, "(1, 2) (2, 0.1) (6, 7) \n(-3, -1) (-2, 0.1) (3, 0.1) ")

  V(real(m1), "1 2 6 \n-3 -2 3 ")
  V(m1.real(), "1 2 6 \n-3 -2 3 ")

  V(imag(m1), "2 0.1 7 \n-1 0.1 0.1 ")
  V(m1.imag(), "2 0.1 7 \n-1 0.1 0.1 ")


  V(herm(m1), "(1, -2) (-3, 1) \n(2, -0.1) (-2, -0.1) \n(6, -7) (3, -0.1) ")
  V(m1.herm(), "(1, -2) (-3, 1) \n(2, -0.1) (-2, -0.1) \n(6, -7) (3, -0.1) ")

  V(trans(m1), "(1, 2) (-3, -1) \n(2, 0.1) (-2, 0.1) \n(6, 7) (3, 0.1) ")
  V(m1.trans(), "(1, 2) (-3, -1) \n(2, 0.1) (-2, 0.1) \n(6, 7) (3, 0.1) ")

  V(conj(m1), "(1, -2) (2, -0.1) (6, -7) \n(-3, 1) (-2, -0.1) (3, -0.1) ")
  V(m1.conj(), "(1, -2) (2, -0.1) (6, -7) \n(-3, 1) (-2, -0.1) (3, -0.1) ")

  VME(m1, herm(herm(m1)))

  VME(herm(m1), conj(trans(m1)))

  VME(conj(m1), herm(trans(m1)))

  VME(trans(m1), herm(conj(m1)))

  VME(m1, trans(herm(conj(m1))))

  double a[] = {1., 2., 3., 4.};
  dmatrix m3(a, 2, 2);

  V(m3, "1 3 \n2 4 ")

  V(conj(m3), "1 3 \n2 4 ")

  V(trans(m3), "1 2 \n3 4 ")

  V(herm(m3), "1 2 \n3 4 ")

END_TEST(MatrixHermTransConj)


BEGIN_TEST(MatrixBlockReplace, "MatrixBlockReplace")
  double a[] = {1., 2., 3., 4., 5., 6.};
  double b[] = {20, 30, 40, 50};
  dmatrix m1(a, 2, 3);
  dmatrix m2(b, 2, 2);

  V(m1.block(0,1,0,2), "1 3 ")
  V(m1.block(0,1,0,1), "1 ")

  m1.replace(0,1, m2);
  V(m1, "1 20 40 \n2 30 50 ")

END_TEST(MatrixBlockReplace)


BEGIN_TEST(MatrixGetRowCol, "MatrixGetRowCol")
  double a[] = {1., 2., 3., 4., 5., 6.};
  dmatrix m1(a, 2, 3);

  V(m1.col(1), "3 \n4 ")

  V(m1.row(1), "2 4 6 ")

END_TEST(MatrixGetRowCol)


#ifdef USE_MKL
BEGIN_TEST(MatrixSolver, "MatrixSolver")

  double a[] = {1., 2., 3., 4., 5., 6., 8, 20, 33};
  dmatrix m1(a, 3, 3);
  
  cvector v;
  dmatrix vl, vr;

  m1.eigen(v, vl, vr);
  V(v, "{(37.6431, 0), (-0.0563885, 0), (1.41334, 0)}")
  V(vl, "-0.090437 -0.451975 -0.0714581 \n-0.190746 0.807446 -0.834439 \n-0.977465 -0.379144 0.546447 ")
  V(vr, "-0.236039 -0.989503 -0.891627 \n-0.51814 0.128295 -0.421717 \n-0.82208 0.0665147 0.164791 ")

  m1.leigen(v, vl);
  V(v, "{(37.6431, 0), (-0.0563885, 0), (1.41334, 0)}")
  V(vl, "-0.090437 -0.451975 -0.0714581 \n-0.190746 0.807446 -0.834439 \n-0.977465 -0.379144 0.546447 ")

  m1.reigen(v, vr);
  V(v, "{(37.6431, 0), (-0.0563885, 0), (1.41334, 0)}")
  V(vr, "-0.236039 -0.989503 -0.891627 \n-0.51814 0.128295 -0.421717 \n-0.82208 0.0665147 0.164791 ")

/* complex */
  _CD b[] = {{1., 1}, {2., 3.2}, {3.,4}, {4.,1}, {5.,1}, {6.,22}, {8,23}, {20,1},{33,2}};
  cmatrix m2(b, 3, 3);

  cvector v1;
  cmatrix vl1, vr1;

  m2.eigen(v1, vl1, vr1);
  V(v1, "{(38.6085, 15.7802), (4.08146, -2.43421), (-3.68992, -9.34598)}")
  V(vl1, "(0.0903046, -0.103983) (-0.0911386, -0.47054) (-0.114558, 0.0915536) \n(0.338112, -0.412648) (0.852899, 0) (0.862411, 0) \n(0.834527, 0) (-0.199099, 0.0566482) (-0.439391, -0.204148) ")
  V(vr1, "(0.342437, 0.337502) (0.979951, 0) (0.82562, 0) \n(0.40302, -0.101577) (0.0410844, 0.0457229) (0.00844035, -0.458344) \n(0.772066, 0) (-0.100078, -0.160943) (-0.319397, 0.078653) ")

  m2.leigen(v1, vl1);
  V(v1, "{(38.6085, 15.7802), (4.08146, -2.43421), (-3.68992, -9.34598)}")
  V(vl1, "(0.0903046, -0.103983) (-0.0911386, -0.47054) (-0.114558, 0.0915536) \n(0.338112, -0.412648) (0.852899, 0) (0.862411, 0) \n(0.834527, 0) (-0.199099, 0.0566482) (-0.439391, -0.204148) ")

  m2.reigen(v1, vr1);
  V(v1, "{(38.6085, 15.7802), (4.08146, -2.43421), (-3.68992, -9.34598)}")
  V(vr1, "(0.342437, 0.337502) (0.979951, 0) (0.82562, 0) \n(0.40302, -0.101577) (0.0410844, 0.0457229) (0.00844035, -0.458344) \n(0.772066, 0) (-0.100078, -0.160943) (-0.319397, 0.078653) ")

 END_TEST(MatrixSolver)
#endif


/**
 * Tests for identity matrix
*/
BEGIN_TEST(IdentityMatrix, "IdentityMatrix")
  dmatrix m1 = didmatrix(3);
  V(m1, "1 0 0 \n0 1 0 \n0 0 1 ")

  cidmatrix m2(2);
  V(m2, "(1, 0) (0, 0) \n(0, 0) (1, 0) ")

END_TEST(IdentityMatrix)


/**
 * Tests for vector 
*/
BEGIN_TEST(VectorTests, "VectorTests")
  double cv1[] = {1, 2, 3, 4};
  double cv2[] = {1, 3, 2, 4};
  
  dvector v1 (cv1, 4);
  dvector v2 (cv2, 4);
  
  double cm1[] = {1, 2, 3, 4, 5, 6};
  dmatrix m1(cm1, 2, 3); 

  V(dvector(m1), "{1, 2, 3, 4, 5, 6}")

  dvector v4;
  v4 = m1;
  V(dvector(v4), "{1, 2, 3, 4, 5, 6}")

  V(dvector(v1), "{1, 2, 3, 4}")
 
  dvector v3;
  v3 = v1;
  V(dvector(v3), "{1, 2, 3, 4}")

  dvector vv(4);
  vv = 1;
  VDE(vv.size(), 4);
  V(vv, "{1, 1, 1, 1}");

  v4.resize(4);
  V(v4, "{1, 2, 3, 4}")

  v4.resize(7,1);
  V(v4, "{1, 2, 3, 4}")

  VDE(v4.size(), 4)

  dvector v5;
  v5 = v1+v2;
  V(v5, "{2, 5, 5, 8}")

  v5 = v1;
  v5 += v2;
  V(v5, "{2, 5, 5, 8}")

  v5 = v1;
  v5 -= v2;
  V(v5, "{0, -1, 1, 0}")

  v5 = v1;
  v5 = v5*2;
  V(v5, "{2, 4, 6, 8}")

  v5 = v1;
  v5 *= 2;
  V(v5, "{2, 4, 6, 8}")

  v5 = v1+v1;
  v5 = v5/2;
  V(v5, "{1, 2, 3, 4}")

  v5 = v1+v1;
  v5 /= 2;
  V(v5, "{1, 2, 3, 4}")

  VDE(max(v5), 4)
  VDE(v5.max(), 4)

  VDE(norm(v5)*norm(v5), 30)
  VDE(v5.norm()*v5.norm(), 30)

  VDE(sum(v5), 10)
  VDE(v5.sum(), 10)

  VDE(v1*v1, 30)
  
  dmatrix m2 = didmatrix(4);
  V(v1*m2, "{1, 2, 3, 4}")
  V(m2*v1, "{1, 2, 3, 4}")

  _CD c1[]= { {1.,2.}, {-3., -1.}, {2.0, 0.1}, {-2., 0.1}, {6., 7.}, {3.0, 0.1}};
  cvector v6(c1, 6);

  V(v6, "{(1, 2), (-3, -1), (2, 0.1), (-2, 0.1), (6, 7), (3, 0.1)}")

  V(real(v6), "{1, -3, 2, -2, 6, 3}")
  V(v6.real(), "{1, -3, 2, -2, 6, 3}")

  V(imag(v6), "{2, -1, 0.1, 0.1, 7, 0.1}")
  V(v6.imag(), "{2, -1, 0.1, 0.1, 7, 0.1}")

  V(trans(v6), "{(1, 2), (-3, -1), (2, 0.1), (-2, 0.1), (6, 7), (3, 0.1)}")

  V(conj(v6), "{(1, -2), (-3, 1), (2, -0.1), (-2, -0.1), (6, -7), (3, -0.1)}")

  V(herm(v6), "{(1, -2), (-3, 1), (2, -0.1), (-2, -0.1), (6, -7), (3, -0.1)}")

  VME(v6, conj(herm(v6)))

  VME(herm(v6), conj(trans(v6)))

  double a[] = {1., 2., 3., 4.};
  dvector m3(a, 4);

  V(m3, "{1, 2, 3, 4}")

  V(conj(m3), "{1, 2, 3, 4}")

  V(trans(m3), "{1, 2, 3, 4}")

  V(herm(m3), "{1, 2, 3, 4}")

  V(v6.block(1,4), "{(-3, -1), (2, 0.1), (-2, 0.1)}")

  cvector v7;
  v7 = v6.trans();
  V(trans(v7), "{(1, 2), (-3, -1), (2, 0.1), (-2, 0.1), (6, 7), (3, 0.1)}")

  v7 = v6.conj();
  V(v7, "{(1, -2), (-3, 1), (2, -0.1), (-2, -0.1), (6, -7), (3, -0.1)}")

  v7 = v6.herm();
  V(v7, "{(1, -2), (-3, 1), (2, -0.1), (-2, -0.1), (6, -7), (3, -0.1)}")
  
  VDE(norm(v6)*norm(v6), 117.03)
  VDE(v6.norm()*v6.norm(), 117.03)

  V(dvector(v1.stl()), "{1, 2, 3, 4}")

END_TEST(VectorTests)

RUN_TEST

  return 0;
}
