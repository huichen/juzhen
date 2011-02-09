ifeq (${MLCPP_BLASLIB}, mkl) 
CFLAGS_MLCPP = -I../ -DUSE_MKL -lmkl_rt -liomp5 -lpthread -DNDEBUG -DMLCPP -O3 -msse2
CC_MLCPP = g++ 
else
ifeq (${MLCPP_BLASLIB}, gotoblas2) 
CFLAGS_MLCPP = -I../ -lgoto2 -DNDEBUG -DMLCPP -O3 -msse2
CC_MLCPP = g++ 
else
CFLAGS_MLCPP = -I../ -lblas -DNDEBUG -DMLCPP -O3 -msse2
CC_MLCPP = g++ 
endif
endif

CFLAGS_EIGEN = -I../../eigen -DNDEBUG -O3 -msse2
CC_EIGEN = g++ 

all:
	time -p $(CC_MLCPP) bench.cpp -o mlcpp $(CFLAGS_MLCPP) 
	time -p $(CC_EIGEN) bench.cpp -o eigen $(CFLAGS_EIGEN) 

run:
	./mlcpp
	./eigen

clean:
	rm *.plt *.txt *.png *.tar mlcpp eigen -rf

plt:
	gnuplot *.plt
