/*
+---------------------------------------------------------------------------+
|  Matrix Library for C++ (mlcpp)                                           |
+---------------------------------------------------------------------------+
|                                                                           |
|  Copyright 2011 Hui Chen                                                  |
|                                                                           |
|  Licensed under the Apache License, Version 2.0 (the "License");          |
|  you may not use this file except in compliance with the License.         |
|  You may obtain a copy of the License at                                  |
|                                                                           |
|      http://www.apache.org/licenses/LICENSE-2.0                           |
|                                                                           |
|  Unless required by applicable law or agreed to in writing, software      |
|  distributed under the License is distributed on an "AS IS" BASIS,        |
|  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. |
|  See the License for the specific language governing permissions and      |
|  limitations under the License.                                           |
|                                                                           |
+---------------------------------------------------------------------------+
*/

#include "mlcpp_test.h"

#include <iostream>
#include <vector>
#include <sstream>

#include <mlcpp.h>

using mlcpp::CS;
using mlcpp::CD;

using mlcpp::smatrix;
using mlcpp::dmatrix;
using mlcpp::cmatrix;
using mlcpp::zmatrix;

using mlcpp::identity_smatrix;
using mlcpp::identity_dmatrix;
using mlcpp::identity_cmatrix;
using mlcpp::identity_zmatrix;

using mlcpp::svector;
using mlcpp::dvector;
using mlcpp::cvector;
using mlcpp::zvector;

int main() {

/**
 * Tests for matrix
*/
BEGIN_TEST(MatrixAssignment, "MatrixAssignment")
  double a[] = {1., 2., 3., 4.};
  double b[] = {1., 2., 3., 4., 5., 6.};

  smatrix m1(a, 2, 2);
  V(m1, "1 3 \n2 4 ")

  smatrix m2 = m1;
  V(m2, "1 3 \n2 4 ")

  m2 = dmatrix(a, 2, 2);
  V(m2, "1 3 \n2 4 ")

  dmatrix mm;
  mm = m2;
  V(mm, "1 3 \n2 4 ")

  m2(0,1) = 99.;
  V(m2, "1 99 \n2 4 ")

  m1 = m2;
  V(m1, "1 99 \n2 4 ")

  smatrix m3(b, 3, 2);
  V(m3, "1 4 \n2 5 \n3 6 ")

  smatrix m4(3,3);
  m4 = 3;
  V(m4, "3 3 3 \n3 3 3 \n3 3 3 ")

  cmatrix m5 = m3;
  V(m5, "(1, 0) (4, 0) \n(2, 0) (5, 0) \n(3, 0) (6, 0) ")

  cmatrix m6;
  m6 = m3;
  V(m6, "(1, 0) (4, 0) \n(2, 0) (5, 0) \n(3, 0) (6, 0) ")

  cmatrix m7(m6.raw_ptr(),2,2);
  V(m7, "(1, 0) (3, 0) \n(2, 0) (4, 0) ")

  smatrix m8(4,4);
  m8 = 3;
  V(m8, "3 3 3 3 \n3 3 3 3 \n3 3 3 3 \n3 3 3 3 ")

  m8 = m8;
  V(m8, "3 3 3 3 \n3 3 3 3 \n3 3 3 3 \n3 3 3 3 ")

  smatrix m9(3,3);
  m9 = 4;
  m9.Clear();
  V(m9, "0 0 0 \n0 0 0 \n0 0 0 ")

END_TEST(MatrixAssignment)

BEGIN_TEST(MatrixElementType, "MatrixElementType")

{
  zmatrix m1(3,3);
  m1 = 1;
  zmatrix m2(3,3);
  m2 = 1;
  zmatrix m3(3,3);
  m3 = 1;
  
  m1 = m2*m3+m1/2;
  V(m1, "(3.5, 0) (3.5, 0) (3.5, 0) \n(3.5, 0) (3.5, 0) (3.5, 0) \n(3.5, 0) (3.5, 0) (3.5, 0) ")
  V(Conjugate(m1)*Adjoint(m1), "(36.75, 0) (36.75, 0) (36.75, 0) \n(36.75, 0) (36.75, 0) (36.75, 0) \n(36.75, 0) (36.75, 0) (36.75, 0) ")
}

{
  cmatrix m1(3,3);
  m1 = 1;
  cmatrix m2(3,3);
  m2 = 1;
  cmatrix m3(3,3);
  m3 = 1;
  
  m1 = m2*m3 + m1/2.0;
  V(m1, "(3.5, 0) (3.5, 0) (3.5, 0) \n(3.5, 0) (3.5, 0) (3.5, 0) \n(3.5, 0) (3.5, 0) (3.5, 0) ")
  V(Conjugate(m1)*Adjoint(m1), "(36.75, 0) (36.75, 0) (36.75, 0) \n(36.75, 0) (36.75, 0) (36.75, 0) \n(36.75, 0) (36.75, 0) (36.75, 0) ")

}

{
  dmatrix m1(3,3);
  m1 = 1;
  dmatrix m2(3,3);
  m2 = 1;
  dmatrix m3(3,3);
  m3 = 1;
  
  m1 = m2*m3+m1/2.0;
  V(m1, "3.5 3.5 3.5 \n3.5 3.5 3.5 \n3.5 3.5 3.5 ")
  V(Conjugate(m1)*Adjoint(m1), "36.75 36.75 36.75 \n36.75 36.75 36.75 \n36.75 36.75 36.75 ")
}

{
  smatrix m1(3,3);
  m1 = 1;
  smatrix m2(3,3);
  m2 = 1;
  smatrix m3(3,3);
  m3 = 1;
  
  m1 = m2*m3+m1/2.0;
  V(m1, "3.5 3.5 3.5 \n3.5 3.5 3.5 \n3.5 3.5 3.5 ")
  V(Conjugate(m1)*Adjoint(m1), "36.75 36.75 36.75 \n36.75 36.75 36.75 \n36.75 36.75 36.75 ")

}


END_TEST(MatrixElementType)

BEGIN_TEST(MatrixCompare, "MatrixCompare")

  double a[] = {1., 2., 3., 4.};
  double b[] = {1., 2., 3., 4., 5., 6.};

  smatrix m1(a, 2, 2);
  smatrix m2(a, 2, 3);
  
  VB(m1==m1)

  VB(m1!=m2)

  smatrix m3 = m1;
  VB(m1==m3)

  m3(1, 1) = 5;
  VB(m1!=m3)

  cmatrix m4 = m1;
/* 
  It's not comparable between cmatrix and smatrix
  VB(m4!=m1);
*/

  cmatrix m5;
  m5 = m4;
  VB(m4==m5)

  m4(1,1) = CS{3,3};
  VB(m4!=m5)

END_TEST(MatrixCompare)

BEGIN_TEST(MatrixPrivateData, "MatrixPrivateData")
  double a[] = {1., 2., 3., 4., 5., 6.};
  smatrix m1(a, 2, 3);

  VDE(m1.num_col(), 3)

  VDE(m1.num_row(), 2)

END_TEST(MatrixPrivateData)

BEGIN_TEST(MatrixResize, "MatrixResize")
  double a[] = {1., 2., 3., 4., 5., 6.};
  smatrix m1(a, 2, 3);
 
  m1.Resize(3,2); 
  V(m1, "1 4 \n2 5 \n3 6 ")
 
  m1.Resize(3,1); 
  V(m1, "1 \n2 \n3 ")

  m1.Resize(3,2); 
  V(m1, "1 4 \n2 5 \n3 6 ")

  m1.Resize(10,10); 
  m1.Resize(3,2); 
  V(m1, "1 4 \n2 5 \n3 6 ")

END_TEST(MatrixResize)

BEGIN_TEST(MatrixIndex, "MatrixIndex")
  double a[] = {1., 2., 3., 4., 5., 6.};
  smatrix m1(a, 2, 3);

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
  smatrix m1(a, 2, 3);

  /**
   * smatrix m2 = m1+m1 will give wrong results. Never use it. 
   */
  smatrix m2;
  m2 = m1+m1;
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
  smatrix m3(b, 2, 2);

  smatrix m4;
  m4 = m3*m2;
  V(m4, "7 15 23 \n10 22 34 ")

  m4 = m3;
  m4 *= m2;
  V(m4, "7 15 23 \n10 22 34 ")

  CS c1[]= { {1.,2.}, {-3., -1.}, {2.0, 0.1}, {-2., 0.1}, {6., 7.}, {3.0, 0.1}};
  cmatrix m5(c1, 2,3);
  V(m5, "(1, 2) (2, 0.1) (6, 7) \n(-3, -1) (-2, 0.1) (3, 0.1) ")

  cmatrix m6;
  m6 = m5 + m5;
  V(m6, "(2, 4) (4, 0.2) (12, 14) \n(-6, -2) (-4, 0.2) (6, 0.2) ")

END_TEST(MatrixArithmatic)

BEGIN_TEST(MatrixHermTransConj, "MatrixHermTransConj")

  CS c1[]= { {1.,2.}, {-3., -1.}, {2.0, 0.1}, {-2., 0.1}, {6., 7.}, {3.0, 0.1}};
  cmatrix m1(c1, 2,3);
  V(m1, "(1, 2) (2, 0.1) (6, 7) \n(-3, -1) (-2, 0.1) (3, 0.1) ")

  V(Real(Real(m1)), "1 2 6 \n-3 -2 3 ")
  V(Imag(Real(m1)), "0 0 0 \n0 0 0 ")

  V(Real(Imag(m1)), "2 0.1 7 \n-1 0.1 0.1 ")


  V(Adjoint(m1), "(1, -2) (-3, 1) \n(2, -0.1) (-2, -0.1) \n(6, -7) (3, -0.1) ")
  V(m1.Adjoint(), "(1, -2) (-3, 1) \n(2, -0.1) (-2, -0.1) \n(6, -7) (3, -0.1) ")

  V(Transpose(m1), "(1, 2) (-3, -1) \n(2, 0.1) (-2, 0.1) \n(6, 7) (3, 0.1) ")
  V(m1.Transpose(), "(1, 2) (-3, -1) \n(2, 0.1) (-2, 0.1) \n(6, 7) (3, 0.1) ")

  V(Conjugate(m1), "(1, -2) (2, -0.1) (6, -7) \n(-3, 1) (-2, -0.1) (3, -0.1) ")
  V(m1.Conjugate(), "(1, -2) (2, -0.1) (6, -7) \n(-3, 1) (-2, -0.1) (3, -0.1) ")

  VME(m1, Adjoint(Adjoint(m1)))

  VME(Adjoint(m1), Conjugate(Transpose(m1)))

  VME(Conjugate(m1), Adjoint(Transpose(m1)))

  VME(Transpose(m1), Adjoint(Conjugate(m1)))

  VME(m1, Transpose(Adjoint(Conjugate(m1))))

  cmatrix m2 = m1;
  V(m2.SwapCol(1,2), "(1, 2) (6, 7) (2, 0.1) \n(-3, -1) (3, 0.1) (-2, 0.1) " )
  m2 = m1;
  m2.SwapCol(1,2);
  V(m2, "(1, 2) (6, 7) (2, 0.1) \n(-3, -1) (3, 0.1) (-2, 0.1) " )

  m2 = m1;
  m2.SwapRow(0,1);
  V(m2, "(-3, -1) (-2, 0.1) (3, 0.1) \n(1, 2) (2, 0.1) (6, 7) " )
  m2 = m1;
  V(m2.SwapRow(0,1), "(-3, -1) (-2, 0.1) (3, 0.1) \n(1, 2) (2, 0.1) (6, 7) " )

  double a[] = {1., 2., 3., 4.};
  smatrix m3(a, 2, 2);

  V(m3, "1 3 \n2 4 ")

  V(Conjugate(m3), "1 3 \n2 4 ")

  V(Transpose(m3), "1 2 \n3 4 ")

  V(Adjoint(m3), "1 2 \n3 4 ")

END_TEST(MatrixHermTransConj)


BEGIN_TEST(MatrixBlockReplace, "MatrixBlockReplace")
  double a[] = {1., 2., 3., 4., 5., 6.};
  double b[] = {20, 30, 40, 50};
  smatrix m1(a, 2, 3);
  smatrix m2(b, 2, 2);

  V(m1.Block(0,1,0,2), "1 3 ")
  V(m1.Block(0,1,0,1), "1 ")
  V(m1.Block(0,1,0,1)*3, "3 ")

  m1.Replace(0,1, m2);
  V(m1, "1 20 40 \n2 30 50 ")

  V((2*m1).Replace(0,1, 2*m2)/2, "1 20 40 \n2 30 50 ")

END_TEST(MatrixBlockReplace)


BEGIN_TEST(MatrixGetRowCol, "MatrixGetRowCol")
  double a[] = {1., 2., 3., 4., 5., 6.};
  smatrix m1(a, 2, 3);

  V(m1.GetColumn(1), "3 \n4 ")

  V(m1.GetRow(1), "2 4 6 ")

END_TEST(MatrixGetRowCol)


#ifdef USE_MKL
BEGIN_TEST(MatrixSolver, "MatrixSolver")

  double a[] = {1., 2., 3., 4., 5., 6., 8, 20, 33};
  dmatrix m1(a, 3, 3);
  
  zvector v;
  dmatrix vl, vr;

  m1.EigenSolver(v, vl, vr);
  V(v, "{(37.6431, 0), (-0.0563885, 0), (1.41334, 0)}")
  V(vl, "-0.090437 -0.451975 -0.0714581 \n-0.190746 0.807446 -0.834439 \n-0.977465 -0.379144 0.546447 ")
  V(vr, "-0.236039 -0.989503 -0.891627 \n-0.51814 0.128295 -0.421717 \n-0.82208 0.0665147 0.164791 ")

  (3*m1).EigenSolver(v, vl, vr);
  V(v/3, "{(37.6431, 0), (-0.0563885, 0), (1.41334, 0)}")
  V(vl, "-0.090437 -0.451975 -0.0714581 \n-0.190746 0.807446 -0.834439 \n-0.977465 -0.379144 0.546447 ")
  V(vr, "-0.236039 -0.989503 -0.891627 \n-0.51814 0.128295 -0.421717 \n-0.82208 0.0665147 0.164791 ")

  m1.LeftEigenSolver(v, vl);
  V(v, "{(37.6431, 0), (-0.0563885, 0), (1.41334, 0)}")
  V(vl, "-0.090437 -0.451975 -0.0714581 \n-0.190746 0.807446 -0.834439 \n-0.977465 -0.379144 0.546447 ")

  m1.RightEigenSolver(v, vr);
  V(v, "{(37.6431, 0), (-0.0563885, 0), (1.41334, 0)}")
  V(vr, "-0.236039 -0.989503 -0.891627 \n-0.51814 0.128295 -0.421717 \n-0.82208 0.0665147 0.164791 ")

/* complex */
  CD b[] = {{1., 1}, {2., 3.2}, {3.,4}, {4.,1}, {5.,1}, {6.,22}, {8,23}, {20,1},{33,2}};
  zmatrix m2(b, 3, 3);

  zvector v1;
  zmatrix vl1, vr1;

  m2.EigenSolver(v1, vl1, vr1);
  V(v1, "{(38.6085, 15.7802), (4.08146, -2.43421), (-3.68992, -9.34598)}")
  V(vl1, "(0.0903046, -0.103983) (-0.0911386, -0.47054) (-0.114558, 0.0915536) \n(0.338112, -0.412648) (0.852899, 0) (0.862411, 0) \n(0.834527, 0) (-0.199099, 0.0566482) (-0.439391, -0.204148) ")
  V(vr1, "(0.342437, 0.337502) (0.979951, 0) (0.82562, 0) \n(0.40302, -0.101577) (0.0410844, 0.0457229) (0.00844035, -0.458344) \n(0.772066, 0) (-0.100078, -0.160943) (-0.319397, 0.078653) ")

  m2.LeftEigenSolver(v1, vl1);
  V(v1, "{(38.6085, 15.7802), (4.08146, -2.43421), (-3.68992, -9.34598)}")
  V(vl1, "(0.0903046, -0.103983) (-0.0911386, -0.47054) (-0.114558, 0.0915536) \n(0.338112, -0.412648) (0.852899, 0) (0.862411, 0) \n(0.834527, 0) (-0.199099, 0.0566482) (-0.439391, -0.204148) ")

  m2.RightEigenSolver(v1, vr1);
  V(v1, "{(38.6085, 15.7802), (4.08146, -2.43421), (-3.68992, -9.34598)}")
  V(vr1, "(0.342437, 0.337502) (0.979951, 0) (0.82562, 0) \n(0.40302, -0.101577) (0.0410844, 0.0457229) (0.00844035, -0.458344) \n(0.772066, 0) (-0.100078, -0.160943) (-0.319397, 0.078653) ")

 END_TEST(MatrixSolver)
#endif


/**
 * Tests for identity matrix
*/
BEGIN_TEST(IdentityMatrix, "IdentityMatrix")
  smatrix m1 = identity_smatrix(3);
  V(m1, "1 0 0 \n0 1 0 \n0 0 1 ")

  identity_cmatrix m2(2);
  V(m2, "(1, 0) (0, 0) \n(0, 0) (1, 0) ")

END_TEST(IdentityMatrix)


/**
 * Tests for vector 
*/
BEGIN_TEST(VectorTests, "VectorTests")
  double cv1[] = {1, 2, 3, 4};
  double cv2[] = {1, 3, 2, 4};
  
  svector v1 (cv1, 4);
  svector v2 (cv2, 4);
  
  double cm1[] = {1, 2, 3, 4, 5, 6};
  smatrix m1(cm1, 2, 3); 

  V(svector(m1), "{1, 2, 3, 4, 5, 6}")

  svector v4;
  v4 = m1;
  V(svector(v4), "{1, 2, 3, 4, 5, 6}")

  V(svector(v1), "{1, 2, 3, 4}")
 
  svector v3;
  v3 = v1;
  V(svector(v3), "{1, 2, 3, 4}")

  svector vv(4);
  vv = 1;
  VDE(vv.size(), 4);
  V(vv, "{1, 1, 1, 1}");

  v4.Resize(4);
  V(v4, "{1, 2, 3, 4}")

  v4.Resize(4,1);
  V(v4, "{1, 2, 3, 4}")

  VDE(v4.size(), 4)

  svector v5;
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

  VDE(Max(v5), 4)
  VDE(v5.Max(), 4)

  VDE(norm(v5)*norm(v5), 30.0)

  VDE(Sum(v5), 10)
  VDE(v5.Sum(), 10)

  VDE(v1*v1, 30)
  
  smatrix m2 = identity_smatrix(4);
  V(v1*m2, "{1, 2, 3, 4}")
  V(m2*v1, "{1, 2, 3, 4}")

  m2(3,3) = 2;
  v1*=m2;
  V(v1, "{1, 2, 3, 8}")

  CS c1[]= { {1.,2.}, {-3., -1.}, {2.0, 0.1}, {-2., 0.1}, {6., 7.}, {3.0, 0.1}};
  cvector v6(c1, 6);

  V(v6, "{(1, 2), (-3, -1), (2, 0.1), (-2, 0.1), (6, 7), (3, 0.1)}")

  V(Real(v6), "{1, -3, 2, -2, 6, 3}")

  V(Imag(v6), "{2, -1, 0.1, 0.1, 7, 0.1}")

  V(Transpose(v6), "{(1, 2), (-3, -1), (2, 0.1), (-2, 0.1), (6, 7), (3, 0.1)}")

  V(Conjugate(v6), "{(1, -2), (-3, 1), (2, -0.1), (-2, -0.1), (6, -7), (3, -0.1)}")

  V(Adjoint(v6), "{(1, -2), (-3, 1), (2, -0.1), (-2, -0.1), (6, -7), (3, -0.1)}")

  VME(v6, Conjugate(Adjoint(v6)))

  VME(Adjoint(v6), Conjugate(Transpose(v6)))

  double a[] = {1., 2., 3., 4.};
  svector m3(a, 4);

  V(m3, "{1, 2, 3, 4}")

  V(Conjugate(m3), "{1, 2, 3, 4}")

  V(Transpose(m3), "{1, 2, 3, 4}")

  V(Adjoint(m3), "{1, 2, 3, 4}")

  V(v6.Block(1,4), "{(-3, -1), (2, 0.1), (-2, 0.1)}")
 
  cvector v7;
  v7 = v6.Transpose();
  V(Transpose(v7), "{(1, 2), (-3, -1), (2, 0.1), (-2, 0.1), (6, 7), (3, 0.1)}")

  v7 = v6.Conjugate();
  V(v7, "{(1, -2), (-3, 1), (2, -0.1), (-2, -0.1), (6, -7), (3, -0.1)}")

  v7 = v6.Adjoint();
  V(v7, "{(1, -2), (-3, 1), (2, -0.1), (-2, -0.1), (6, -7), (3, -0.1)}")
  
  VSE(norm(v6)*norm(v6), 117.03)

  V(svector(v1.STL()), "{1, 2, 3, 8}")

  svector v8 =  m3;
  V(v8.Swap(1,2), "{1, 3, 2, 4}")

  v8 = m3;
  v8.Swap(1,2);
  V(v8, "{1, 3, 2, 4}")

END_TEST(VectorTests)

RUN_TEST

  return 0;
}
