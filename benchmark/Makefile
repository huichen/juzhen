ifeq (${MLCPP_BLASLIB}, mkl) 
CFLAGS_MLCPP = -I../src -DUSE_MKL -lmkl_rt -liomp5 -lpthread -DNDEBUG -DMLCPP -O3 -msse2
CC_MLCPP = g++ 
else
ifeq (${MLCPP_BLASLIB}, gotoblas2) 
CFLAGS_MLCPP = -DUSE_BLAS -I../src -lgoto2 -fopenmp -DNDEBUG -DMLCPP -O3 -msse2
CC_MLCPP = g++ 
else
ifeq (${MLCPP_BLASLIB}, blas) 
CFLAGS_MLCPP = -DUSE_BLAS -I../src -lblas -DNDEBUG -DMLCPP -O3 -msse2
CC_MLCPP = g++ 
else
CFLAGS_MLCPP = -I../src -DNDEBUG -DMLCPP -O3 -msse2 -fopenmp
CC_MLCPP = g++ 
endif
endif
endif

CFLAGS_EIGEN = -I../../eigen -DNDEBUG -O3 -msse2 -fopenmp
CC_EIGEN = g++ 

all: benchmark.cc benchmark.h
	time -p $(CC_MLCPP) benchmark.cc -o juzhen $(CFLAGS_MLCPP)
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

