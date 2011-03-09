const double a1[] = {1, -2, -3, 4};

dmatrix m1(a1, 2, 2);
zmatrix m2(a1, 2, 2);

V(Abs(m1), "1 3 \n2 4 ")
V(Abs2(m1), "1 9 \n4 16 ")

V(Abs(m2), "1 3 \n2 4 ")
V(Abs2(m2), "1 9 \n4 16 ")

const double a2[] = {1, -2, -3, 4, 7, 2, 8, 9, 7};
dmatrix m3(a2, 3, 3);
zmatrix m4(a2, 3, 3);
#ifdef USE_MKL
VDE(Det(m3), 115)
VDE(Det(m4).real, 115)
#endif
