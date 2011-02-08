ifeq (${MLCPP_BLASLIB}, mkl) 
CFLAGS_MLCPP = -I../ -DUSE_MKL -lmkl_rt -liomp5 -lpthread -DNDEBUG -DMLCPP -O3
CC_MLCPP = g++ 
else
CFLAGS_MLCPP = -I../ -lblas 
CC_MLCPP = g++ 
endif

CFLAGS_EIGEN = -I../../eigen -DNDEBUG -O3 
CC_EIGEN = g++ 

all:
	time -p $(CC_MLCPP) bench.cpp -o mlcpp $(CFLAGS_MLCPP) 
	time -p $(CC_EIGEN) bench.cpp -o eigen $(CFLAGS_EIGEN) 
	./mlcpp
	./eigen

clean:
	rm *.plt *.txt mlcpp eigen -rf

plt:
	gnuplot *.plt
