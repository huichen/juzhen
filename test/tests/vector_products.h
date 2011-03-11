dvector v1(3);
v1 << 1, 2, 3;
dvector v2(3);
v2 << 2, 1, 7;

V(CrossProduct(v1, v2), "{11, -1, -3}")

V(OuterProduct(v1, v2), "2 1 7 \n4 2 14 \n6 3 21 ")
