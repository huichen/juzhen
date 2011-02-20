ifeq (${MLCPP_BLASLIB}, mkl) 
CFLAGS = -DUSE_MKL -I../include -lmkl_rt -liomp5 -lpthread -std=c++0x
CC = g++ 
else
ifeq (${MLCPP_BLASLIB}, gotoblas2) 
CFLAGS = -I../include -lgoto2 -pthread -std=c++0x 
CC = g++ 
else
CFLAGS = -I../include -lblas -std=c++0x 
CC = g++ 
endif
endif

all: schrodinger

schrodinger: schrodinger.cc
	$(CC) schrodinger.cc -o schrodinger $(CFLAGS) -O3 -DNDEBUG 
