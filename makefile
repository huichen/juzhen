all:
	icc test.cpp -o test -lmkl_rt -g
	./test

clean:
	rm test
