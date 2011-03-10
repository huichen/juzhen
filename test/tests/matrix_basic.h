CS c1[]= { {1., 2.}, {-3., -1.}, {2.0, 0.1}, {-2., 0.1}, {6., 7.}, {3.0, 0.1}};
cmatrix m1(c1, 2, 3);
VDE(m1.num_col(), 3)
VDE(m1.num_row(), 2)

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
V(m2.SwapCol(1, 2), "(1, 2) (6, 7) (2, 0.1) \n(-3, -1) (3, 0.1) (-2, 0.1) " )
m2 = m1;
m2.SwapCol(1, 2);
V(m2, "(1, 2) (6, 7) (2, 0.1) \n(-3, -1) (3, 0.1) (-2, 0.1) " )

m2 = m1;
m2.SwapRow(0, 1);
V(m2, "(-3, -1) (-2, 0.1) (3, 0.1) \n(1, 2) (2, 0.1) (6, 7) " )
m2 = m1;
V(m2.SwapRow(0, 1), "(-3, -1) (-2, 0.1) (3, 0.1) \n(1, 2) (2, 0.1) (6, 7) " )

double a[] = {1., 2., 3., 4.};
smatrix m3(a, 2, 2);

V(m3, "1 3 \n2 4 ")

V(Conjugate(m3), "1 3 \n2 4 ")

V(Transpose(m3), "1 2 \n3 4 ")

V(Adjoint(m3), "1 2 \n3 4 ")

const double a1[] = {1, -2, -3, 4};

dmatrix m4(a1, 2, 2);
zmatrix m5(a1, 2, 2);

V(Abs(m4), "1 3 \n2 4 ")
V(Abs2(m4), "1 9 \n4 16 ")

V(Abs(m5), "1 3 \n2 4 ")
V(Abs2(m5), "1 9 \n4 16 ")

const double a2[] = {1, -2, -3, 4, 7, 2, 8, 9, 7};
dmatrix m6(a2, 3, 3);
zmatrix m7(a2, 3, 3);
#ifdef USE_MKL
VDE(Det(m6), 115)
VDE(Det(m7).real, 115)
#endif
