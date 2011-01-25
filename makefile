all:
	icc test.cpp -o test -lmkl_rt
	./test
