#include <iostream>
#include <mklpp.h>
#include <unittest.hpp>
#include <vector>
#include <sstream>

using namespace mklpp;

int main() {

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

  m4(1,1) = (MKLCD){3,3};
  VB(m4!=m5)

END_TEST(MatrixCompare)

BEGIN_TEST(MatrixPrivateData, "MatrixPrivateData")
  double a[] = {1., 2., 3., 4., 5., 6.};
  dmatrix m1(a, 2, 3);

  VDE(m1.nCol(), 3)

  VDE(m1.nRow(), 2)

  VDE(m1.getDataSize(), 6)

  VDE(m1.getTranspose(), CblasNoTrans)

  m1.trans();
  VDE(m1.getTranspose(), CblasTrans)

  m1.trans();
  VDE(m1.getTranspose(), CblasNoTrans)

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


BEGIN_TEST(MatrixHermTransConj, "MatrixHermTransConj")

  MKLCD c1[]= { {1.,2.}, {-3., -1.}, {2.0, 0.}, {-2., .0}, {6., 7.}, {3.0, 0.}};
  cmatrix m1(c1, 2,3);
  V(m1, "(1, 2) (2, 0) (6, 7) \n(-3, -1) (-2, 0) (3, 0) ")

  V(herm(m1), "(1, -2) (-3, 1) \n(2, 0) (-2, 0) \n(6, -7) (3, 0) ")

  V(trans(m1), "(1, 2) (-3, -1) \n(2, 0) (-2, 0) \n(6, 7) (3, 0) ")

  V(conj(m1), "(1, -2) (2, 0) (6, -7) \n(-3, 1) (-2, 0) (3, 0) ")

  VME(m1, herm(herm(m1)))

  VME(herm(m1), conj(trans(m1)))

  VME(conj(m1), herm(trans(m1)))

  VME(trans(m1), herm(conj(m1)))

  VME(m1, trans(herm(conj(m1))))

  cmatrix m2;
  m2 = herm(m1);
  V(m2, "(1, -2) (-3, 1) \n(2, 0) (-2, 0) \n(6, -7) (3, 0) ")

  m2 = trans(m1);
  V(m2, "(1, 2) (-3, -1) \n(2, 0) (-2, 0) \n(6, 7) (3, 0) ")

  m2 = conj(m1);
  V(m2, "(1, -2) (2, 0) (6, -7) \n(-3, 1) (-2, 0) (3, 0) ")


END_TEST(MatrixHermTransConj)

BEGIN_TEST(IdentityMatrix, "IdentityMatrix")
  dmatrix m1 = didmatrix(3);
  V(m1, "1 0 0 \n0 1 0 \n0 0 1 ")

  cidmatrix m2(2);
  V(m2, "(1, 0) (0, 0) \n(0, 0) (1, 0) ")

END_TEST(IdentityMatrix)


BEGIN_TEST(VectorOperation, "VectorOperation")
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

  v4.resize(4);
  V(v4, "{1, 2, 3, 4}")

  v4.resize(7,1);
  V(v4, "{1, 2, 3, 4}")

  VDE(v4.size(), 4)

  VDE(v1*v1, 30)
  
  V(dvector(v1.stl()), "{1, 2, 3, 4}")


END_TEST(VectorOperation)

  UnitTest MyTest;
  MyTest.Run(); 

  return 0;
}
