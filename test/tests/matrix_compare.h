smatrix m1(2, 2);
m1 << 1, 3, 2, 4;
smatrix m2(2, 3);
m2 << 1, 3, 5,
      2, 4, 6;

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

m4(1, 1) = CS(3, 3);
VB(m4!=m5)
