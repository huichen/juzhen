%pythoncode {

def inverse(a):
    return a.inverse()

def linear_solver(a, b):
    return a.linear_solver(b)

def eigen(a):
    e = CMatrix()
    if type(a).__name__ == "Matrix":
      vl = Matrix()
      vr = Matrix()
    else:
      vl = CMatrix()
      vr = CMatrix()
    a.eigen(e, vl, vr)
    return (e, vl, vr)

def right_eigen(a):
    e = CMatrix()
    if type(a).__name__ == "Matrix":
      vr = Matrix()
    else:
      vr = CMatrix()
    a.right_eigen(e, vr)
    return (e, vr)

def left_eigen(a):
    e = CMatrix()
    if type(a).__name__ == "Matrix":
      vr = Matrix()
    else:
      vr = CMatrix()
    a.left_eigen(e, vr)
    return (e, vr)

}
