all:
	icc mklcpp.cpp -o mklcpp -lmkl_rt
	./mklcpp
