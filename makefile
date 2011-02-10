ifeq (${MLCPP_BLASLIB}, mkl) 
CFLAGS = -DUSE_MKL -I./ -lmkl_rt -liomp5 -lpthread -std=c++0x
CC = g++ 
else
ifeq (${MLCPP_BLASLIB}, gotoblas2) 
CFLAGS = -I./ -lgoto2 -pthread -std=c++0x 
CC = g++ 
else
CFLAGS = -I./ -lblas -std=c++0x 
CC = g++ 
endif
endif

test:
	$(CC) unittest.cpp -o unittest $(CFLAGS) -g
	./unittest

sample:
	$(CC) example.cpp -o example $(CFLAGS) -O3 -DNDEBUG 
	./example

clean:
	rm -rf unittest example html latex

doc:
	doxygen resources/doxygen.config
