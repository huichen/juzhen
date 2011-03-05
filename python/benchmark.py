#! /usr/bin/env python

from numpy import *
from mlcpp import CMatrix

import time

min_rank = 10
max_rank = 500
step_rank = 10

def worker_numpy(n):
  t0 = time.time()
  a = ones((n, n), dtype = complex128)
  b = ones((n, n), dtype = complex128)
  for i in range(10):
    c = dot(a, b)
  t1 = time.time()
  return t1 - t0

def worker_mlcpp(n):
  t0 = time.time()
  a = CMatrix(n, n)
  a.set(1)
  b = CMatrix(n, n)
  b.set(1)
  for i in range(10):
    c = a * b
  t1 = time.time()
  return t1 - t0

print "Rank\tmlcpp(s)\t\tnumpy(s)"
f = open("benchmark.txt", "w")
for i in range(min_rank, max_rank, step_rank):
  print '%d\t%6f\t%6f' % (i, worker_mlcpp(i), worker_numpy(i))
  f.write('%d\t%6f\t%6f\n' % (i, worker_mlcpp(i), worker_numpy(i)))
f.close()

