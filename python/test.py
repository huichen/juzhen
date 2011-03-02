from mlpy import Complex, Conjugate, abs, abs2, real, imag

a = Complex(3, 4)
b = Complex(-1)
c = Complex(b)

print a
print b 
print c 

print Conjugate(a)
print abs2(a)
print abs(a)
print real(a)
print imag(a)

print a + b
print a + 3

print a - b
print a - 3

print a * b
print a * 3

print a / b
print a / 3
