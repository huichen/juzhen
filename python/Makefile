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
	swig -c++ -python mlpy.i
	$(CC) -c -fpic mlpy_wrap.cxx -I/usr/include/python2.6/ $(CFLAGS)
	$(CC) -shared mlpy_wrap.o -o _mlpy.so $(CFLAGS)
	python test.py > output
	
clean:
	rm *.o mlpy*.pyc *_wrap* -rf

distclean: clean
	rm *.so mlpy.py -rf
