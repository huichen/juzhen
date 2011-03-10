ifeq (${JUZHEN_BLASLIB}, mkl) 
CFLAGS_JUZHEN = -I../src -DUSE_MKL -lmkl_rt -liomp5 -lpthread -DNDEBUG -DJUZHEN -O3 -msse2
CC_JUZHEN = g++ 
else
ifeq (${JUZHEN_BLASLIB}, gotoblas2) 
CFLAGS_JUZHEN = -DUSE_BLAS -I../src -lgoto2 -fopenmp -DNDEBUG -DJUZHEN -O3 -msse2
CC_JUZHEN = g++ 
else
ifeq (${JUZHEN_BLASLIB}, blas) 
CFLAGS_JUZHEN = -DUSE_BLAS -I../src -lblas -DNDEBUG -DJUZHEN -O3 -msse2
CC_JUZHEN = g++ 
else
CFLAGS_JUZHEN = -I../src -DNDEBUG -DJUZHEN -O3 -msse2 -fopenmp
CC_JUZHEN = g++ 
endif
endif
endif

CFLAGS_EIGEN = -I../../eigen -DNDEBUG -O3 -msse2 -fopenmp
CC_EIGEN = g++ 

all: benchmark.cc benchmark.h
	time -p $(CC_JUZHEN) benchmark.cc -o juzhen $(CFLAGS_JUZHEN)
	time -p $(CC_EIGEN) benchmark.cc -o eigen $(CFLAGS_EIGEN)

run_small:
	./juzhen 100 2 1 small_
	./eigen 100 2 1 small_

run_big:
	./juzhen 1000 100 10 big_
	./eigen 1000 100 10 big_

clean:
	rm *.plt *.txt *.png *.tar -rf

plt:
	gnuplot *.plt

