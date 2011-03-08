double c1[] = {1, 2, 3, 4};
double c2[] = {1, 5, 2, 4};

dmatrix m1(c1, 2, 2); 
dmatrix m2(c2, 2, 2); 

VDE(Trace(m1), 5);
V(SchurProduct(m1, m2), "1 6 \n10 16 ")
V(KroneckerProduct(m1, m2), "1 2 3 6 \n5 4 15 12 \n2 4 4 8 \n10 8 20 16 ")
VDE(InnerProduct(m1, m2), 33)
VDE(InnerProduct(m1), 30)

