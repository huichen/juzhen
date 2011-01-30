all:
	icc quantum.cpp -o quantum -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -pthread -liomp5 -O3 
	./quantum

t:
	icc test.cpp -o test -lmkl_rt -liomp5 -lpthread -g
	./test

clean:
	rm test quantum
