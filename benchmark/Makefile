ifeq (${MLCPP_BLASLIB}, mkl) 
CFLAGS_MLCPP = -I../include -DUSE_MKL -lmkl_rt -liomp5 -lpthread -DNDEBUG -DMLCPP -O3 -msse2
CC_MLCPP = g++ 
else
ifeq (${MLCPP_BLASLIB}, gotoblas2) 
CFLAGS_MLCPP = -I../include -lgoto2 -DNDEBUG -DMLCPP -O3 -msse2
CC_MLCPP = g++ 
else
CFLAGS_MLCPP = -I../include -lblas -DNDEBUG -DMLCPP -O3 -msse2
CC_MLCPP = g++ 
endif
endif

CFLAGS_EIGEN = -I../../eigen -DNDEBUG -O3 -msse2
CC_EIGEN = g++ 

all:
	time -p $(CC_MLCPP) benchmark.cc -o mlcpp $(CFLAGS_MLCPP)
	time -p $(CC_EIGEN) benchmark.cc -o eigen $(CFLAGS_EIGEN)

run:
	./mlcpp
	./eigen

clean:
	rm *.plt *.txt *.png *.tar mlcpp eigen -rf

plt:
	gnuplot *.plt

