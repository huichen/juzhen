test:
	g++ unittest.cpp -DUSE_ATLAS -I./ -o unittest -lblas -llapack -g
	./unittest

sample:
	g++ example.cpp -DUSE_ATLAS -I./ -o example -lblas -llapack 
	./example

clean:
	rm -rf unittest example 
