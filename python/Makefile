ifeq (${JUZHEN_BLASLIB}, mkl) 
CFLAGS = -DUSE_MKL -I../src -lmkl_rt -liomp5 -lpthread -lmkl_intel_thread -lmkl_core -std=c++0x -O3 -msse2 -DNDEBUG
CC = g++ 
else
ifeq (${JUZHEN_BLASLIB}, gotoblas2) 
CFLAGS = -I../src -lgoto2 -pthread -std=c++0x -O3 -msse2 -DNDEBUG
CC = g++ 
else
ifeq (${JUZHEN_BLASLIB}, blas) 
CFLAGS = -I../src -lblas -std=c++0x -O3 -msse2 -DNDEBUG
CC = g++ 
else
CFLAGS = -I../src -std=c++0x -O3 -msse2 -DNDEBUG -fopenmp
CC = g++ 
endif
endif
endif

all: juzhen_wrap.cc util.i
	$(CC) -c -fpic juzhen_wrap.cc -I/usr/include/python2.6/ $(CFLAGS)
	$(CC) -shared juzhen_wrap.o -o _juzhen.so $(CFLAGS)
	
wrap:
	swig -c++ -python  -o juzhen_wrap.cc juzhen.i

util:
	./extract_util.py

clean:
	rm *.o juzhen*.pyc -rf

distclean: clean
	rm *.so juzhen.py *.txt *.png -rf

bench: _juzhen.so juzhen.py benchmark.py
	./benchmark.py

plt: benchmark.txt benchmark.plt
	gnuplot benchmark.plt	
	
