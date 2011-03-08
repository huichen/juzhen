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

m4(1, 1) = CS{3, 3};
VB(m4!=m5)
