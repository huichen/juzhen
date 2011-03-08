double c1[] = {1., 2., 3.};
double c2[] = {2., 1., 7.};
dvector v1(c1, 3);
dvector v2(c2, 3);

V(CrossProduct(v1, v2), "{11, -1, -3}")

V(OuterProduct(v1, v2), "2 1 7 \n4 2 14 \n6 3 21 ")
