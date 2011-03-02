#! /usr/bin/env python

from mlpy import Complex, Matrix, CMatrix

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
a.Clear()
for i in range(0, 3):
  for j in range(0, 3):
    a.Set(i, j, i*j+1)
print "a = "
print a
print

print "a.Get(1, 1) = "
print a.Get(1, 1)
print

b = a.Copy();
b.Set(1, 1, 3)
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

print "a.GetCol(2) = "
print a.GetCol(2)
print

print "a.GetRow(1) = "
print a.GetRow(1)
print

print "a.Block(1, 3, 0, 2) = "
print a.Block(1, 3, 0, 2)
print

c = Matrix(2, 2)
c.Clear()
print "a.Replace(1, 0, c) = "
print a.Replace(1, 0, c)
print

print "a.SwapCol(0, 1) = "
print a.SwapCol(0, 1)
print

print "a.SwapRow(0, 1) = "
print a.SwapRow(0, 1)
print

print "a.Resize(2, 2); a = "
a.Resize(2, 2)
print a
print

print "========= Test CMatrix ========="
a = CMatrix(3, 3)
a.Clear()
for i in range(0, 3):
  for j in range(0, 3):
    a.Set(i, j, i*j+1)
print "a = "
print a
print

print "a.Get(1, 1) = "
print a.Get(1, 1)
print

b = a.Copy();
b.Set(1, 1, 3)
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

