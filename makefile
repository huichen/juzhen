all:
	icc example.cpp -o example -lmkl_rt -pthread -liomp5 -O3 
	./example

test:
	icc unittest.cpp -o unittest -lmkl_rt -liomp5 -lpthread -g
	./unittest

clean:
	rm -rf unittest example 
