const double a1[] = {1, -2, -3, 4};

dmatrix m1(a1, 2, 2);
zmatrix m2(a1, 2, 2);

V(Abs(m1), "1 3 \n2 4 ")
V(Abs2(m1), "1 9 \n4 16 ")

V(Abs(m2), "1 3 \n2 4 ")
V(Abs2(m2), "1 9 \n4 16 ")
