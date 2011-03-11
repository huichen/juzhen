dvector v1(4);
v1 << 1, 2, 3, 4;
dvector v2(4);
v2 << 1, 3, 2, 4;

VDE(Max(v1), 4)
VDE(Min(v1), 1)
VDE(Sum(v1), 10)
VDE(Prod(v1), 24)
VDE(Mean(v1), 2.5)

VDE(Norm(v2)*Norm(v2), 30.0)
VDE(NormSquare(v2), 30.0)
VDE(NormOne(v2), 10.0)
VDE(NormInfinity(v2), 4)

V(Sort(v2), "{1, 2, 3, 4}")

VDE(Cov(v1, v2), 1)
VDE(Var(v1), 1.25)
VDE(StdDev(v1), 1.11803)
VDE(CorrCoeff(v1, v2), 0.8) 
