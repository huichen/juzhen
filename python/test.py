#! /usr/bin/env python

from mlpy import Complex, Matrix, CMatrix, Vector, CVector

print "========= Test Complex ========="
a = Complex(3, 4)
b = Complex(-1)
c = Complex(b)

print "a = ", a
print "b = ", b 
print "c = ",c 

print "a.conj() = ", a.conj()
print "a.abs() = ", a.abs() 
print "a.real = ", a.real
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
a = Matrix(3, 3)
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

print "a.get_col(2) = "
print a.get_col(2)
print

print "a.get_row(1) = "
print a.get_row(1)
print

print "a.block(1, 3, 0, 2) = "
print a.block(1, 3, 0, 2)
print

c = Matrix(2, 2)
c.clear()
print "a.replace(1, 0, c) = "
print a.replace(1, 0, c)
print

print "a.swap_col(0, 1) = "
print a.swap_col(0, 1)
print

print "a.swap_row(0, 1) = "
print a.swap_row(0, 1)
print

print "a.resize(2, 2); a = "
a.resize(2, 2)
print a
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
a = Vector(6)
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

print "a.sum() = "
print a.sum()
print 

print "a.max() = "
print a.max()
print 

print "a.norm() = "
print a.norm()
print 

print "b.sort() = "
print b.sort()
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

print "a.sum() = "
print a.sum()
print 

print "a.norm() = "
print a.norm()
print 
