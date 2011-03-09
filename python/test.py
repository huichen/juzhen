#! /usr/bin/env python

from mlcpp import Complex, Matrix, CMatrix, Identity, CIdentity, Vector, CVector
from mlcpp import inverse 
from mlcpp import linear_solver
from mlcpp import eigen, left_eigen, right_eigen

print "========= Test Complex ========="
# a = 3 + 4i
a = Complex(3, 4)

# a = -1 + 0i 
b = Complex(-1)

# This is equivalent to c = b
c = Complex(b)

print "a = ", a
print "b = ", b 
print "c = ", c 

# Conjugate of a
print "a.conj() = ", a.conj()

# Absolute value of a
print "a.abs() = ", a.abs() 

# Real part of a
print "a.real = ", a.real

# Imaginary part of a
print "a.imag = ", a.imag 

print "a + b = ", a + b
print "a + 3 = ", a + 3

print "a - b = ", a - b
print "a - 3 = ", a - 3

print "a * b = ", a * b
print "a * 3 = ", a * 3

print "a / b = ", a / b
print "a / 3 = ", a / 3

print "========= Test Matrix ========="
# Define a matrix of 3 rows (1st parameter) and 3 columns (2nd parameter)
# Row and column indices start from zero
# 
#     x    x    x    (row 0)
#     x    x    x    (row 1)
#     x    x    x    (row 2)
#
# Col 0    1    2
# 
a = Matrix(3, 3)

# Number of rows and columns
print "a.num_row() = ", a.num_row()
print "a.num_col() = ", a.num_col()

# Set zero to a's elements
a.clear()

# Set 3 to all elements
a.set(3)
print "a.set(3); a = "
print a
print

# 7 x 7 identity matrix
print "Identity(7) = "
print Identity(7)
print

for i in range(0, 3):
  for j in range(0, 3):
    a.set(i, j, i*j+1) # Set a's element at i-th row & j-th row to be i*j+1
print "a = "
print a
print

# Get a's (1,2) element
print "a.get(1, 2) = "
print a.get(1, 2)
print

# Get a copy of a and save to b
# b = a doesn't copy the values
b = a.copy();
b.set(1, 1, 3)
print "b = "
print b
print

print "a + b = "
print a + b
print 

print "a - b = "
print a - b
print 

print "a * b = "
print a * b
print 

print "a * 3 = "
print a * 3 
print 

print "a / 3 = "
print a / 3 
print 

# Get column 2
print "a.get_col(2) = "
print a.get_col(2)
print

# Get row 1
print "a.get_row(1) = "
print a.get_row(1)
print

# Get a block between row 1 (include) to row 3 (not include) and column 0 (include) to column 2 (not include). 
print "a.block(1, 3, 0, 2) = "
print a.block(1, 3, 0, 2)
print

c = Matrix(2, 2)
c.clear()

# Replace the sub matrix with upper-left corner at (1, 0) with matrix c
print "a.replace(1, 0, c) = "
print a.replace(1, 0, c)
print

# Swap column 0 and column 1 
print "a.swap_col(0, 1) = "
print a.swap_col(0, 1)
print

# Swap row 0 and row 1 
print "a.swap_row(0, 1) = "
print a.swap_row(0, 1)
print

# Resize the dimension of a
print "a.resize(2, 2); a = "
a.resize(2, 2)
print a
print

# Get the real part of a matrix
print "a.real() = "
print a.real()
print

# Get the imaginary part of a matrix
print "a.imag() = "
print a.imag()
print

# Get the transpose of a matrix
print "a.trans() = "
print a.trans()
print

# Get the conjugate of a matrix
print "a.conj() = "
print a.conj()
print

# Get the adjoint of a matrix
print "a.adj() = "
print a.adj()
print

print "========= Test CMatrix ========="
a = CMatrix(3, 3)
a.clear()
for i in range(0, 3):
  for j in range(0, 3):
    a.set(i, j, i*j+1)
print "a = "
print a
print

print "a.get(1, 1) = "
print a.get(1, 1)
print

b = a.copy();
b.set(1, 1, 3)
print "b = "
print b
print

print "a + b = "
print a + b
print 

print "a - b = "
print a - b
print 

print "a * b = "
print a * b
print 

print "a * 3 = "
print a * 3 
print 

print "a / 3 = "
print a / 3 
print 

print "a.real() = "
print a.real() 
print 

print "a.imag() = "
print a.imag() 
print 

print "========= Test Vector ========="
# Define a vector of length 6
a = Vector(6)

# Set zero to all elements of a
a.clear()

# Set 3 to all elements
a.set(3)
print "a = "
print a
print

# Set elements
for i in range(0, 6):
    a.set(i, i)
print "a = "
print a
print

# Get size of a vector
print "a.size() = "
print a.size()
print

# Get element 1 
print "a.get(1) = "
print a.get(1)
print

# Copy data of a to vector b
b = a.copy();

# Set element 1 to be 3
b.set(1, 3)
print "b = "
print b
print

print "a + b = "
print a + b
print 

print "a - b = "
print a - b
print 

print "a * b = "
print a * b
print 

print "a * 3 = "
print a * 3 
print 

print "a / 3 = "
print a / 3 
print 

print "a * Identity(6) = "
print a * Identity(6)
print 

print "Identity(6) * a= "
print Identity(6) * a
print 

# Get the real part of a vector
print "a.real() = "
print a.real()
print

# Get the imaginary part of a vector
print "a.imag() = "
print a.imag()
print

print "========= Test CVector ========="
a = CVector(6)
a.clear()
for i in range(0, 6):
    a.set(i, i)
print "a = "
print a
print

print "a.get(1) = "
print a.get(1)
print

b = a.copy();
b.set(1, 3)
print "b = "
print b
print

print "a + b = "
print a + b
print 

print "a - b = "
print a - b
print 

print "a * b = "
print a * b
print 

print "a * 3 = "
print a * 3 
print 

print "a / 3 = "
print a / 3 
print 

print "a.real() = "
print a.real()
print

print "a.imag() = "
print a.imag()
print

print "========= Test Matrix Solver ========="

a = Identity(3)
print "inverse(a) = "
print inverse(a)
print

b = Matrix(3, 3)
b.set(10)

# Solve A * X = B
print "linear_solver(a, b) = "
print linear_solver(a, b)
print

a.set(1, 2, 3)
a.set(0, 1, 7)
a.set(2, 0, 7)

# Solve eigen problem:
#   e is a column vector contains all eigen values
#   vl is a square matrix that contains all left eigen vectors 
#   vr is a square matrix that contains all right eigen vectors 
# Similarly,
# e, vl = left_eigen(a) solves left eigen vectors only
# e, vr = right_eigen(a) solves right eigen vectors only
e, vl, vr = eigen(a);
print "e, vl, vr = eigen(a); e = "
print e 
print
print "vl = "
print vl 
print
print "vr = "
print vr 
print

