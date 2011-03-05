%rename(conj) Matrix::Conjugate;
%rename(adj) Matrix::Adjoint;
%rename(trans) Matrix::Transpose;
%rename(resize) Matrix::Resize;
%rename(clear) Matrix::Clear;
%rename(block) Matrix::Block;
%rename(replace) Matrix::Replace;
%rename(swap_col) Matrix::SwapCol;
%rename(swap_row) Matrix::SwapRow;
%rename(get_col) Matrix::GetCol;
%rename(get_row) Matrix::GetRow;

%rename(conj) CMatrix::Conjugate;
%rename(adj) CMatrix::Adjoint;
%rename(trans) CMatrix::Transpose;
%rename(resize) CMatrix::Resize;
%rename(clear) CMatrix::Clear;
%rename(block) CMatrix::Block;
%rename(replace) CMatrix::Replace;
%rename(swap_col) CMatrix::SwapCol;
%rename(swap_row) CMatrix::SwapRow;
%rename(get_col) CMatrix::GetCol;
%rename(get_row) CMatrix::GetRow;

class Matrix {
 public:
  Matrix();
  Matrix(int, int);
  Matrix(const Matrix &);
  ~Matrix();

  bool operator==(const Matrix &);
  bool operator!=(const Matrix &);

  int num_row();
  int num_col();

  void Resize(int, int);
  void Clear();

  Matrix Conjugate();
  Matrix Transpose();
  Matrix Adjoint();

  Matrix GetCol(int);
  Matrix GetRow(int);

  Matrix Block(int, int, int, int);
  Matrix &Replace(int, int, const Matrix &);

  Matrix &SwapCol(int, int);
  Matrix &SwapRow(int, int);

  %extend {
    char *__str__() {
      static char temp[1024];
      snprintf(temp, 1024, "%s", OutputToString(*$self).c_str());
      return temp;
    }

    double &__getitem__(unsigned int i) {
      return $self->operator[](i);
    }

    Matrix copy() {
      return *$self;
    }

    Matrix &set(int r, int c, double d) {
      $self->operator()(r, c) = d;
      return *$self;
    }

    Matrix &set(int r, double d) {
      $self->operator()(r) = d;
      return *$self;
    }

    Matrix &set(double d) {
      $self->Set(d);
      return *$self;
    }

    double get(int r, int c) {
      return $self->operator()(r, c);
    }

    double get(int r) {
      return $self->operator()(r);
    }

    Matrix __add__(const Matrix &a) {
      return (*$self) + a;
    }

    Matrix __sub__(const Matrix &a) {
      return (*$self) - a;
    }

    Matrix __mul__(const Matrix &a) {
      return (*$self) * a;
    }

    Vector __mul__(const Vector &a) {
      return (*$self) * a;
    }

    Matrix __mul__(double a) {
      return (*$self) * a;
    }

    Matrix __div__(double a) {
      return (*$self) / a;
    }

    Matrix real() {
      return Real(*$self);
    }

    Matrix imag() {
      return Imag(*$self);
    }

    Matrix inverse() {
      return Inverse(*$self);
    }

    Matrix linear_solver(const Matrix &b) {
      return LinearSolver(*$self, b);
    }

    void eigen(CMatrix &e, Matrix &vl, Matrix &vr) {
      EigenSolver(*$self, e, vl, vr);
    }

    void left_eigen(CMatrix &e, Matrix &vl) {
      LeftEigenSolver(*$self, e, vl);
    }

    void right_eigen(CMatrix &e, Matrix &vr) {
      RightEigenSolver(*$self, e, vr);
    }
  }
}; 

class CMatrix {
 public:
  CMatrix();
  CMatrix(int, int);
  CMatrix(const CMatrix &);
  ~CMatrix();

  bool operator==(const CMatrix &);
  bool operator!=(const CMatrix &);

  int num_row();
  int num_col();

  void Resize(int, int);
  void Clear();

  CMatrix Conjugate();
  CMatrix Transpose();
  CMatrix Adjoint();

  CMatrix GetCol(int);
  CMatrix GetRow(int);

  CMatrix Block(int, int, int, int);
  CMatrix &Replace(int, int, const CMatrix &);

  CMatrix &SwapCol(int, int);
  CMatrix &SwapRow(int, int);

  %extend {
    char *__str__() {
      static char temp[1024];
      snprintf(temp, 1024, "%s", OutputToString(*$self).c_str());
      return temp;
    }

    Complex &__getitem__(unsigned int i) {
      return $self->operator[](i);
    }

    CMatrix copy() {
      return *$self;
    }

    CMatrix &set(int r, int c, Complex &d) {
      $self->operator()(r, c) = d;
      return *$self;
    }

    CMatrix &set(int r, int c, double d) {
      $self->operator()(r, c) = d;
      return *$self;
    }

    CMatrix &set(int r, Complex &d) {
      $self->operator()(r) = d;
      return *$self;
    }

    CMatrix &set(int r, double d) {
      $self->operator()(r) = d;
      return *$self;
    }

    CMatrix &set(Complex &d) {
      $self->Set(d);
      return *$self;
    }

    CMatrix &set(double d) {
      $self->Set(d);
      return *$self;
    }

    Complex get(int r, int c) {
      return $self->operator()(r, c);
    }

    Complex get(int r) {
      return $self->operator()(r);
    }

    CMatrix __add__(const CMatrix &a) {
      return (*$self) + a;
    }

    CMatrix __sub__(const CMatrix &a) {
      return (*$self) - a;
    }

    CMatrix __mul__(const CMatrix &a) {
      return (*$self) * a;
    }

    CMatrix __mul__(Complex &a) {
      return (*$self) * a;
    }

    CVector __mul__(const CVector &a) {
      return (*$self) * a;
    }

    CMatrix __mul__(double a) {
      return (*$self) * a;
    }

    CMatrix __div__(Complex &a) {
      return (*$self) / a;
    }

    CMatrix __div__(double a) {
      return (*$self) / a;
    }

    Matrix real() {
      return Real(*$self);
    }

    Matrix imag() {
      return Imag(*$self);
    }
    
    CMatrix inverse() {
      return Inverse(*$self);
    }

    CMatrix linear_solver(const CMatrix &b) {
      return LinearSolver(*$self, b);
    }

    void eigen(CMatrix &e,CMatrix &vl,CMatrix &vr) {
      EigenSolver(*$self, e, vl, vr);
    }

    void left_eigen(CMatrix &e,CMatrix &vl) {
      LeftEigenSolver(*$self, e, vl);
    }

    void right_eigen(CMatrix &e,CMatrix &vr) {
      RightEigenSolver(*$self, e, vr);
    }
}
}; 

Matrix Identity(int);

CMatrix CIdentity(int);
