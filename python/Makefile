ifeq (${MLCPP_BLASLIB}, mkl) 
CFLAGS = -DUSE_MKL -I../include -lmkl_rt -liomp5 -lpthread -std=c++0x -O3
CC = g++ 
else
ifeq (${MLCPP_BLASLIB}, gotoblas2) 
CFLAGS = -I../include -lgoto2 -pthread -std=c++0x -O3
CC = g++ 
else
CFLAGS = -I../include -lblas -std=c++0x -O3
CC = g++ 
endif
endif

all:
	swig -c++ -python mlpy.i
	$(CC) -c -fpic mlpy_wrap.cxx -I/usr/include/python2.6/ $(CFLAGS) 
	$(CC) -shared mlpy_wrap.o -o _mlpy.so $(CFLAGS) 
	python test.py
	
clean:
	rm *.o *.so mlpy*.py* *_wrap* -rf
