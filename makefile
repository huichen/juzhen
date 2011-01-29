all:
	icc quantum.cpp -o quantum -lmkl_rt -liomp5 -lpthread -O3
	./quantum

test:
	icc test.cpp -o test -lmkl_rt -liomp5 -lpthread -g
	./test

clean:
	rm test quantum
