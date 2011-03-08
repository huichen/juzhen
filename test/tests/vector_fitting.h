double c1[] = {18, 14, 9, 10, 5, 22, 14, 12};
double c2[] = {39, 9, 9, 7, 8, 35, 36, 22};

dvector v1(c1, 8);
dvector v2(c2, 8);

VDE(LeastSquaresMethod(v1, v2), -1.64968)
