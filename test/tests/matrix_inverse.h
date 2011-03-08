#ifdef USE_MKL
// float type
float a[] = {1., 2., 3., 4., 5., 6., 8, 20, 33};
smatrix m1(a, 3, 3);
V(Inverse(m1), "-15 28 -13.3333 \n2 -3 1.33333 \n1 -2 1 ");

// double type
dmatrix m2(a, 3, 3);
V(Inverse(m2), "-15 28 -13.3333 \n2 -3 1.33333 \n1 -2 1 ");

// double type
cmatrix m3(a, 3, 3);
V(Inverse(m3), "(-15, -0) (28, 0) (-13.3333, 0) \n(2, 0) (-3, 0) (1.33333, 0) \n(1, 0) (-2, -0) (1, 0) ");

// double type
zmatrix m4(a, 3, 3);
V(Inverse(m4), "(-15, -0) (28, 0) (-13.3333, 0) \n(2, 0) (-3, 0) (1.33333, 0) \n(1, 0) (-2, -0) (1, 0) ");

#endif
