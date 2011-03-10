%rename(conj) Vector::Conjugate;
%rename(adj) Vector::Adjoint;
%rename(trans) Vector::Transpose;
%rename(resize) Vector::Resize;
%rename(clear) Vector::Clear;
%rename(block) Vector::Block;
%rename(replace) Vector::Replace;
%rename(swap) Vector::Swap;

%rename(conj) CVector::Conjugate;
%rename(adj) CVector::Adjoint;
%rename(trans) CVector::Transpose;
%rename(resize) CVector::Resize;
%rename(clear) CVector::Clear;
%rename(block) CVector::Block;
%rename(replace) CVector::Replace;
%rename(swap) CVector::Swap;

class Vector {
 public:
  Vector();
  Vector(int);
  Vector(const Vector &);
  Vector(const Matrix &);
  ~Vector();

  bool operator==(const Vector &);
  bool operator!=(const Vector &);

  int size();

  void Resize(int);
  void Resize(int, int);
  void Clear();

  Vector Conjugate();
  Vector Transpose();
  Vector Adjoint();

  Vector &Sort(); 

  Vector Block(int, int);
  Vector &Replace(int, const Vector &);

  Vector &Swap(int, int);

  %extend {
    char *__str__() {
      static char temp[1024];
      snprintf(temp, 1024, "%s", String(*$self).c_str());
      return temp;
    }

    double &__getitem__(unsigned int i) {
      return $self->operator[](i);
    }

    Vector copy() {
      return *$self;
    }

    Vector &set(int r, double d) {
      $self->operator()(r) = d;
      return *$self;
    }

    Vector &set(double d) {
      $self->Set(d);
      return *$self;
    }

    double get(int r) {
      return $self->operator()(r);
    }

    Vector __add__(const Vector &a) {
      return (*$self) + a;
    }

    Vector __sub__(const Vector &a) {
      return (*$self) - a;
    }

    Vector __mul__(const Matrix &a) {
      return (*$self) * a;
    }

    double __mul__(const Vector &a) {
      return (*$self) * a;
    }

    Vector __mul__(double a) {
      return (*$self) * a;
    }

    Vector __div__(double a) {
      return (*$self) / a;
    }

    Vector real() {
      return Real(*$self);
    }

    Vector imag() {
      return Imag(*$self);
    }

    Vector abs() {
      return Abs(*$self);
    }

    Vector abs2() {
      return Abs2(*$self);
    }
  }
}; 

class CVector {
 public:
  CVector();
  CVector(int);
  CVector(const CVector &);
  CVector(const CMatrix &);
  ~CVector();

  bool operator==(const CVector &);
  bool operator!=(const CVector &);

  int size();

  void Resize(int);
  void Resize(int, int);
  void Clear();

  CVector Conjugate();
  CVector Transpose();
  CVector Adjoint();

  CVector Block(int, int);
  CVector &Replace(int, const CVector &);

  CVector &Swap(int, int);

  %extend {
    char *__str__() {
      static char temp[1024];
      snprintf(temp, 1024, "%s", String(*$self).c_str());
      return temp;
    }

    Complex &__getitem__(unsigned int i) {
      return $self->operator[](i);
    }

    CVector copy() {
      return *$self;
    }

    CVector &set(int r, double d) {
      $self->operator()(r) = d;
      return *$self;
    }

    CVector &set(int r, Complex d) {
      $self->operator()(r) = d;
      return *$self;
    }

    CVector &set(double d) {
      $self->Set(d);
      return *$self;
    }

    CVector &set(Complex &d) {
      $self->Set(d);
      return *$self;
    }

    Complex get(int r) {
      return $self->operator()(r);
    }

    CVector __add__(const CVector &a) {
      return (*$self) + a;
    }

    CVector __sub__(const CVector &a) {
      return (*$self) - a;
    }

    CMatrix __mul__(const CMatrix &a) {
      return (*$self) * a;
    }

    Complex __mul__(const CVector &a) {
      return (*$self) * a;
    }

    CVector __mul__(double a) {
      return (*$self) * a;
    }

    CVector __mul__(const Complex &a) {
      return $self->operator*(a);
    }

    CVector __div__(double a) {
      return (*$self) / a;
    }

    CVector __div__(const Complex &a) {
      return (*$self) / a;
    }

    Vector real() {
      return Real(*$self);
    }

    Vector imag() {
      return Imag(*$self);
    }

    Vector abs() {
      return Abs(*$self);
    }

    Vector abs2() {
      return Abs2(*$self);
    }
  }
}; 

