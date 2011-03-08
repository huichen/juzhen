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

m2(0, 1) = 99.;
V(m2, "1 99 \n2 4 ")

m1 = -(-m2);
V(m1, "1 99 \n2 4 ")

smatrix m3(b, 3, 2);
V(m3, "1 4 \n2 5 \n3 6 ")

smatrix m4(3, 3);
m4.Set(3);
V(m4, "3 3 3 \n3 3 3 \n3 3 3 ")

cmatrix m5 = m3;
V(m5, "(1, 0) (4, 0) \n(2, 0) (5, 0) \n(3, 0) (6, 0) ")

cmatrix m6;
m6 = m3;
V(m6, "(1, 0) (4, 0) \n(2, 0) (5, 0) \n(3, 0) (6, 0) ")

cmatrix m7(m6.raw_ptr(), 2, 2);
V(m7, "(1, 0) (3, 0) \n(2, 0) (4, 0) ")

smatrix m8(4, 4);
m8.Set(3);
V(m8, "3 3 3 3 \n3 3 3 3 \n3 3 3 3 \n3 3 3 3 ")

m8 = m8;
V(m8, "3 3 3 3 \n3 3 3 3 \n3 3 3 3 \n3 3 3 3 ")

smatrix m9(3, 3);
m9.Set(4);
m9.Clear();
V(m9, "0 0 0 \n0 0 0 \n0 0 0 ")
