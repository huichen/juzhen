/**
 * This file is automatically generated by a script.
 * DON'T MANUALLY CHANGE IT!
 */


/**
 * Functions extracted from '../src/util/matrix_basic.h'.
 */
%rename(transpose) Transpose;
Matrix Transpose(const Matrix &m);
CMatrix Transpose(const CMatrix &m);

%rename(adjoint) Adjoint;
Matrix Adjoint(const Matrix &m);
CMatrix Adjoint(const CMatrix &m);

%rename(conjugate) Conjugate;
Matrix Conjugate(const Matrix &m);
CMatrix Conjugate(const CMatrix &m);

%rename(det) Det;
double Det(const Matrix &m);
Complex Det(const CMatrix &m);


/**
 * Functions extracted from '../src/util/vector_basic.h'.
 */
%rename(transpose) Transpose;
Vector Transpose(const Vector &v);
CVector Transpose(const CVector &v);

%rename(adjoint) Adjoint;
Vector Adjoint(const Vector &v);
CVector Adjoint(const CVector &v);

%rename(conjugate) Conjugate;
Vector Conjugate(const Vector &v);
CVector Conjugate(const CVector &v);

%rename(join) Join;
Vector Join(const Vector &vector1, const Vector &vector2);
CVector Join(const CVector &vector1, const CVector &vector2);


/**
 * Functions extracted from '../src/util/vector_math.h'.
 */
%rename(outer_product) OuterProduct;
Matrix OuterProduct(const Vector &vector1, const Vector &vector2);
CMatrix OuterProduct(const CVector &vector1, const CVector &vector2);

%rename(cross_product) CrossProduct;
Vector CrossProduct(const Vector &vector1, const Vector &vector2);
CVector CrossProduct(const CVector &vector1, const CVector &vector2);


/**
 * Functions extracted from '../src/util/matrix_statistics.h'.
 */
%rename(max) Max;
double Max(const Matrix &matrix);

%rename(min) Min;
double Min(const Matrix &matrix);

%rename(sum) Sum;
double Sum(const Matrix &matrix);
Complex Sum(const CMatrix &matrix);

%rename(prod) Prod;
double Prod(const Matrix &matrix);
Complex Prod(const CMatrix &matrix);

%rename(mean) Mean;
double Mean(const Matrix &matrix);
Complex Mean(const CMatrix &matrix);

%rename(norm_square) NormSquare;
double NormSquare(const Matrix &m);
Complex NormSquare(const CMatrix &m);

%rename(norm) Norm;
double Norm(const Matrix &m);
double Norm(const CMatrix &m);

%rename(norm_one) NormOne;
double NormOne(const Matrix &m);
Complex NormOne(const CMatrix &m);

%rename(norm_infinity) NormInfinity;
double NormInfinity(const Matrix &m);
Complex NormInfinity(const CMatrix &m);


/**
 * Functions extracted from '../src/util/vector_statistics.h'.
 */
%rename(sum) Sum;
double Sum(const Vector &vector);
Complex Sum(const CVector &vector);

%rename(prod) Prod;
double Prod(const Vector &vector);
Complex Prod(const CVector &vector);

%rename(mean) Mean;
double Mean(const Vector &vector);
Complex Mean(const CVector &vector);

%rename(norm_square) NormSquare;
double NormSquare(const Vector &v);
Complex NormSquare(const CVector &v);

%rename(norm) Norm;
double Norm(const Vector &v);
double Norm(const CVector &v);

%rename(norm_one) NormOne;
double NormOne(const Vector &v);
Complex NormOne(const CVector &v);

%rename(norm_infinity) NormInfinity;
double NormInfinity(const Vector &v);
Complex NormInfinity(const CVector &v);

%rename(max) Max;
double Max(const Vector &vector);

%rename(min) Min;
double Min(const Vector &vector);

%rename(sort) Sort;
Vector Sort(const Vector &v);

%rename(cov) Cov;
double Cov(const Vector &vector1, const Vector &vector2);

%rename(var) Var;
double Var(const Vector &vector);

%rename(std_dev) StdDev;
double StdDev(const Vector &vector);

%rename(corr_coeff) CorrCoeff;
double CorrCoeff(const Vector &vector1, const Vector &vector2);


/**
 * Functions extracted from '../src/util/random.h'.
 */
%rename(randomize) Randomize;
Vector &Randomize(Vector &vector, double scale);
CVector &Randomize(CVector &vector, Complex scale);

%rename(randomize) Randomize;
Vector &Randomize(Vector &vector);
CVector &Randomize(CVector &vector);

%rename(randomize) Randomize;
Matrix& Randomize(Matrix &matrix, double scale);
CMatrix& Randomize(CMatrix &matrix, Complex scale);

%rename(randomize) Randomize;
Matrix& Randomize(Matrix &matrix);
CMatrix& Randomize(CMatrix &matrix);


/**
 * Functions extracted from '../src/util/vector_fitting.h'.
 */
%rename(least_squares_method) LeastSquaresMethod;
double LeastSquaresMethod(const Vector &a, const Vector &b);
Complex LeastSquaresMethod(const CVector &a, const CVector &b);


/**
 * Functions extracted from '../src/util/matrix_map_math.h'.
 */
%rename(cos) Cos;
Matrix Cos(const Matrix &matrix);

%rename(sin) Sin;
Matrix Sin(const Matrix &matrix);

%rename(tan) Tan;
Matrix Tan(const Matrix &matrix);

%rename(acos) Acos;
Matrix Acos(const Matrix &matrix);

%rename(asin) Asin;
Matrix Asin(const Matrix &matrix);

%rename(atan) Atan;
Matrix Atan(const Matrix &matrix);

%rename(cosh) Cosh;
Matrix Cosh(const Matrix &matrix);

%rename(sinh) Sinh;
Matrix Sinh(const Matrix &matrix);

%rename(tanh) Tanh;
Matrix Tanh(const Matrix &matrix);

%rename(exp) Exp;
Matrix Exp(const Matrix &matrix);

%rename(log) Log;
Matrix Log(const Matrix &matrix);

%rename(log10) Log10;
Matrix Log10(const Matrix &matrix);

%rename(sqrt) Sqrt;
Matrix Sqrt(const Matrix &matrix);

%rename(ceil) Ceil;
Matrix Ceil(const Matrix &matrix);

%rename(fabs) Fabs;
Matrix Fabs(const Matrix &matrix);

%rename(floor) Floor;
Matrix Floor(const Matrix &matrix);


/**
 * Functions extracted from '../src/util/matrix_math.h'.
 */
%rename(trace) Trace;
double Trace(const Matrix &matrix);
Complex Trace(const CMatrix &matrix);

%rename(schur_product) SchurProduct;
Matrix SchurProduct(const Matrix &matrix1, const Matrix &matrix2);
CMatrix SchurProduct(const CMatrix &matrix1, const CMatrix &matrix2);

%rename(kronecker_product) KroneckerProduct;
Matrix KroneckerProduct(const Matrix &matrix1, const Matrix &matrix2);
CMatrix KroneckerProduct(const CMatrix &matrix1, const CMatrix &matrix2);

%rename(inner_product) InnerProduct;
double InnerProduct(const Matrix &matrix1, const Matrix &matrix2);

%rename(inner_product) InnerProduct;
double InnerProduct(const Matrix &matrix);

%rename(inner_product) InnerProduct;
double InnerProduct(const CMatrix &matrix);


/**
 * Functions extracted from '../src/util/vector_map_math.h'.
 */
%rename(cos) Cos;
Vector Cos(const Vector &vector);

%rename(sin) Sin;
Vector Sin(const Vector &vector);

%rename(tan) Tan;
Vector Tan(const Vector &vector);

%rename(acos) Acos;
Vector Acos(const Vector &vector);

%rename(asin) Asin;
Vector Asin(const Vector &vector);

%rename(atan) Atan;
Vector Atan(const Vector &vector);

%rename(cosh) Cosh;
Vector Cosh(const Vector &vector);

%rename(sinh) Sinh;
Vector Sinh(const Vector &vector);

%rename(tanh) Tanh;
Vector Tanh(const Vector &vector);

%rename(exp) Exp;
Vector Exp(const Vector &vector);

%rename(log) Log;
Vector Log(const Vector &vector);

%rename(log10) Log10;
Vector Log10(const Vector &vector);

%rename(sqrt) Sqrt;
Vector Sqrt(const Vector &vector);

%rename(ceil) Ceil;
Vector Ceil(const Vector &vector);

%rename(fabs) Fabs;
Vector Fabs(const Vector &vector);

%rename(floor) Floor;
Vector Floor(const Vector &vector);


