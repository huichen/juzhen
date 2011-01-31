test:
	icc unittest.cpp -I./ -o unittest -lmkl_rt -liomp5 -lpthread -g
	./unittest

sample:
	icc example.cpp -I./ -o example -lmkl_rt -pthread -liomp5 -O3 
	./example

clean:
	rm -rf unittest example 
