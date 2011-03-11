dmatrix v1(2, 2);
v1 << 1, 3,
      2, 4;

VDE(Max(v1), 4)
VDE(Min(v1), 1)
VDE(Sum(v1), 10)
VDE(Prod(v1), 24)
VDE(Mean(v1), 2.5)

VDE(Norm(v1)*Norm(v1), 30.0)
VDE(NormSquare(v1), 30.0)
VDE(NormOne(v1), 10.0)
VDE(NormInfinity(v1), 4)
