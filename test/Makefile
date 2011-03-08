ifeq (${MLCPP_BLASLIB}, mkl) 
CFLAGS = -DUSE_MKL -I../source -lmkl_rt -liomp5 -lpthread -std=c++0x
CC = g++ 
else
ifeq (${MLCPP_BLASLIB}, gotoblas2) 
CFLAGS = -DUSE_BLAS -I../source -lgoto2 -pthread -std=c++0x
CC = g++ 
else
ifeq (${MLCPP_BLASLIB}, blas) 
CFLAGS = -DUSE_BLAS -I../source -lblas -std=c++0x 
CC = g++ 
else
CFLAGS = -I../source -std=c++0x -fopenmp
CC = g++ 
endif
endif
endif

test: 
	$(CC) mlcpp_test.cc -o unittest $(CFLAGS) -g
	./unittest

clean:
	rm -rf unittest 

