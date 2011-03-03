#! /usr/bin/env python

from numpy import *
from mlpy import CMatrix

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

def worker_mlpy(n):
  t0 = time.time()
  a = CMatrix(n, n)
  a.set(1)
  b = CMatrix(n, n)
  b.set(1)
  for i in range(10):
    c = a * b
  t1 = time.time()
  return t1 - t0

print "Rank\tmlpy(s)\t\tnumpy(s)"
for i in range(min_rank, max_rank, step_rank):
  print '%d\t%6f\t%6f' % (i, worker_mlpy(i), worker_numpy(i))

