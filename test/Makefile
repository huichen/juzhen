ifeq (${MLCPP_BLASLIB}, mkl) 
CFLAGS = -DUSE_MKL -I../include -lmkl_rt -liomp5 -lpthread -std=c++0x
CC = g++ 
else
ifeq (${MLCPP_BLASLIB}, gotoblas2) 
CFLAGS = -I../include -lgoto2 -pthread -std=c++0x 
CC = g++ 
else
CFLAGS = -I../include -lblas -std=c++0x 
CC = g++ 
endif
endif

test: 
	$(CC) mlcpp_test.cc -o unittest $(CFLAGS) -g
	./unittest

clean:
	rm -rf unittest 

