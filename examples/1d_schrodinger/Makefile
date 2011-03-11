ifeq (${JUZHEN_BLASLIB}, mkl) 
CFLAGS = -DUSE_MKL -I../../src -lmkl_rt -liomp5 -lpthread
CC = g++ 
else
ifeq (${JUZHEN_BLASLIB}, gotoblas2) 
CFLAGS = -DUSE_BLAS -I../../src -lgoto2 -pthread 
CC = g++ 
else
ifeq (${JUZHEN_BLASLIB}, gotoblas2) 
CFLAGS = -DUSE_BLAS -I../../src -lblas 
CC = g++ 
else
CFLAGS = -I../../src 
CC = g++ 
endif
endif
endif

all: schrodinger

schrodinger: schrodinger.cc
	$(CC) schrodinger.cc -o schrodinger $(CFLAGS) -O3 -DNDEBUG 
