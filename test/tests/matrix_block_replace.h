double a[] = {1., 2., 3., 4., 5., 6.};
double b[] = {20, 30, 40, 50};
smatrix m1(a, 2, 3);
smatrix m2(b, 2, 2);

V(m1.Block(0, 1, 0, 2), "1 3 ")
V(m1.Block(0, 1, 0, 1), "1 ")
V(m1.Block(0, 1, 0, 1)*3, "3 ")

m1.Replace(0, 1, m2);
V(m1, "1 20 40 \n2 30 50 ")

V((2*m1).Replace(0, 1, 2*m2)/2, "1 20 40 \n2 30 50 ")
