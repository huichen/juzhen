#ifdef USE_MKL
// float type
float a[] = {1., 2., 3., 4., 5., 6., 8, 20, 33};
float b[] = {1., 3., 3., 7., 5., 8., 8, 20, 3};

smatrix m11(a, 3, 3);
smatrix m12(b, 3, 3);

V(LinearSolver(m11, m12), "29 -71.6667 400 \n-3 9.66667 -40 \n-2 5 -29 ");

// double type
dmatrix m21(a, 3, 3);
dmatrix m22(b, 3, 3);

V(LinearSolver(m21, m22), "29 -71.6667 400 \n-3 9.66667 -40 \n-2 5 -29 ");

// single complex type
cmatrix m31(a, 3, 3);
cmatrix m32(b, 3, 3);

V(LinearSolver(m31, m32), "(29, 0) (-71.6667, 0) (400, 0) \n(-3, 0) (9.66667, 0) (-40, 0) \n(-2, -0) (5, 0) (-29, -0) ");

// double complex type
cmatrix m41(a, 3, 3);
cmatrix m42(b, 3, 3);

V(LinearSolver(m41, m42), "(29, 0) (-71.6667, 0) (400, 0) \n(-3, 0) (9.66667, 0) (-40, 0) \n(-2, -0) (5, 0) (-29, -0) ");

#endif
