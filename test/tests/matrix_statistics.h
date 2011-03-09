double c1[] = {1, 2, 3, 4};

dmatrix v1(c1, 2, 2);

VDE(Max(v1), 4)
VDE(Min(v1), 1)
VDE(Sum(v1), 10)
VDE(Prod(v1), 24)
VDE(Mean(v1), 2.5)

VDE(Norm(v1)*Norm(v1), 30.0)
VDE(NormSquare(v1), 30.0)
VDE(NormOne(v1), 10.0)
VDE(NormInfinity(v1), 4)
