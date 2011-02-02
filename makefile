ifeq (${MLCPP_BLASLIB}, mkl) 
CFLAGS = -DUSE_MKL -I./ -lmkl_rt -liomp5 -lpthread 
else
CFLAGS = -I./ -lblas 
endif

test:
	g++ unittest.cpp -o unittest $(CFLAGS) -g
	./unittest

sample:
	g++ example.cpp -o example $(CFLAGS) 
	./example

clean:
	rm -rf unittest example 
