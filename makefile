ifeq (${MLCPP_BLASLIB}, mkl) 
CFLAGS = -DUSE_MKL -I./ -lmkl_rt -liomp5 -lpthread 
CC = icc
else
CFLAGS = -I./ -lblas 
CC = g++ 
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
	doxygen . 
