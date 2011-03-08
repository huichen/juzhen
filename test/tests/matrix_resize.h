double a[] = {1., 2., 3., 4., 5., 6.};
smatrix m1(a, 2, 3);

m1.Resize(3, 2);
V(m1, "1 4 \n2 5 \n3 6 ")

m1.Resize(3, 1);
V(m1, "1 \n2 \n3 ")

m1.Resize(3, 2);
V(m1, "1 4 \n2 5 \n3 6 ")

m1.Resize(10, 10);
m1.Resize(3, 2);
V(m1, "1 4 \n2 5 \n3 6 ")
