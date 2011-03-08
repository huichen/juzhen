double c1[] = {1, 2, 3, 4};
double c2[] = {1, 3, 2, 4};

dvector v1(c1, 4);
dvector v2(c2, 4);

VDE(Max(v1), 4)
VDE(Min(v1), 1)
VDE(Sum(v1), 10)
VDE(Average(v1), 2.5)

VDE(Norm(v2)*Norm(v2), 30.0)
