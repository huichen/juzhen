ifeq (${JUZHEN_BLASLIB}, mkl) 
CFLAGS = -DUSE_MKL -I../src -lmkl_rt -liomp5 -lpthread -std=c++0x
CC = g++ 
else
ifeq (${JUZHEN_BLASLIB}, gotoblas2) 
CFLAGS = -DUSE_BLAS -I../src -lgoto2 -pthread -std=c++0x
CC = g++ 
else
ifeq (${JUZHEN_BLASLIB}, blas) 
CFLAGS = -DUSE_BLAS -I../src -lblas -std=c++0x 
CC = g++ 
else
CFLAGS = -I../src -std=c++0x -fopenmp
CC = g++ 
endif
endif
endif

test: 
	$(CC) test.cc -o unittest $(CFLAGS) -g

gen:
	./generate_tests.py

clean:
	rm -rf unittest 

