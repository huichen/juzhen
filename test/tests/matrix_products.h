
dmatrix m1(2, 2); 
m1 << 1, 3,
      2, 4;
dmatrix m2(2, 2); 
m2 << 1, 2,
      5, 4;

VDE(Trace(m1), 5);
V(SchurProduct(m1, m2), "1 6 \n10 16 ")
V(KroneckerProduct(m1, m2), "1 2 3 6 \n5 4 15 12 \n2 4 4 8 \n10 8 20 16 ")
VDE(InnerProduct(m1, m2), 33)
VDE(InnerProduct(m1), 30)

