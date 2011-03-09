ifeq (${MLCPP_BLASLIB}, mkl) 
CFLAGS = -DUSE_MKL -I../src -lmkl_rt -liomp5 -lpthread -lmkl_intel_thread -lmkl_core -std=c++0x -O3 -msse2 -DNDEBUG
CC = g++ 
else
ifeq (${MLCPP_BLASLIB}, gotoblas2) 
CFLAGS = -I../src -lgoto2 -pthread -std=c++0x -O3 -msse2 -DNDEBUG
CC = g++ 
else
ifeq (${MLCPP_BLASLIB}, blas) 
CFLAGS = -I../src -lblas -std=c++0x -O3 -msse2 -DNDEBUG
CC = g++ 
else
CFLAGS = -I../src -std=c++0x -O3 -msse2 -DNDEBUG -fopenmp
CC = g++ 
endif
endif
endif

all: mlcpp_wrap.cc mlcpp_util.i
	$(CC) -c -fpic mlcpp_wrap.cc -I/usr/include/python2.6/ $(CFLAGS)
	$(CC) -shared mlcpp_wrap.o -o _mlcpp.so $(CFLAGS)
	python test.py > output
	
wrap:
	swig -c++ -python  -o mlcpp_wrap.cc mlcpp.i

util:
	./extract_util.py

clean:
	rm *.o mlcpp*.pyc -rf

distclean: clean
	rm *.so mlcpp.py *.txt *.png -rf

bench: _mlcpp.so mlcpp.py benchmark.py
	./benchmark.py

plt: benchmark.txt benchmark.plt
	gnuplot benchmark.plt	
	
