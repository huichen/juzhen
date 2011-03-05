ifeq (${MLCPP_BLASLIB}, mkl) 
CFLAGS = -DUSE_MKL -I../include -lmkl_rt -liomp5 -lpthread -lmkl_intel_thread -lmkl_core -std=c++0x -O3 -msse2 -DNDEBUG
CC = g++ 
else
ifeq (${MLCPP_BLASLIB}, gotoblas2) 
CFLAGS = -I../include -lgoto2 -pthread -std=c++0x -O3 -msse2 -DNDEBUG
CC = g++ 
else
CFLAGS = -I../include -lblas -std=c++0x -O3 -msse2 -DNDEBUG
CC = g++ 
endif
endif

all:
	swig -c++ -python mlcpp.i
	$(CC) -c -fpic mlcpp_wrap.cxx -I/usr/include/python2.6/ $(CFLAGS)
	$(CC) -shared mlcpp_wrap.o -o _mlcpp.so $(CFLAGS)
	python test.py > output
	
clean:
	rm *.o mlcpp*.pyc *_wrap* -rf

distclean: clean
	rm *.so mlcpp.py *.txt *.png -rf

bench: _mlcpp.so mlcpp.py benchmark.py
	./benchmark.py

plt: benchmark.txt benchmark.plt
	gnuplot benchmark.plt	
	
