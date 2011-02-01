test:
	g++ unittest.cpp -I./ -o unittest -lblas -llapack -g
	./unittest

sample:
	g++ example.cpp -I./ -o example -lblas -llapack 
	./example

clean:
	rm -rf unittest example 
