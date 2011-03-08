ifeq (${MLCPP_BLASLIB}, mkl) 
CFLAGS = -DUSE_MKL -I../src -lmkl_rt -liomp5 -lpthread -std=c++0x
CC = g++ 
else
ifeq (${MLCPP_BLASLIB}, gotoblas2) 
CFLAGS = -DUSE_BLAS -I../src -lgoto2 -pthread -std=c++0x 
CC = g++ 
else
ifeq (${MLCPP_BLASLIB}, gotoblas2) 
CFLAGS = -DUSE_BLAS -I../src -lblas -std=c++0x 
CC = g++ 
else
CFLAGS = -I../src -std=c++0x 
CC = g++ 
endif
endif
endif

all: schrodinger

schrodinger: schrodinger.cc
	$(CC) schrodinger.cc -o schrodinger $(CFLAGS) -O3 -DNDEBUG 
