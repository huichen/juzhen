#ifndef BENCH_HPP
#define BENCH_HPP

#ifdef MLCPP
#include <mlcpp.hpp>
#else
#include <Eigen/Dense>
#endif

#include <fstream>
#include <string.h>
#include <sstream>
#include <sys/time.h>
#include <iostream>

#define N 400
#define NB 10

static int counter=0;

#ifdef MLCPP 
#define FNAME "mlcpp_"
#else
#define FNAME "eigen_"
#endif

#define BEGIN_BENCH(s) \
{\
std::stringstream out;\
out << counter;\
filename =out.str()+ ".plt";\
myfile.open(filename.c_str());\
myfile << "set term png \nset out \"" + out.str()+".png\"\nset title \"" + s + "\"\nset xlabel \"N\" \nset ylabel \"Seconds\" \nplot 'mlcpp_"+out.str() + ".txt' title 'mlcpp' using 1:2 w l, 'eigen_"+out.str() + ".txt' title 'eigen' using 1:2 w l ";\
myfile.close();\
}\
{\
std::stringstream out;\
out << counter++;\
filename =FNAME+out.str()+ ".txt";\
}\
myfile.open(filename.c_str());\
for (size_t ni=10; ni<=1024; ni+=10) {\
H1.resize(ni,ni);\
H2.resize(ni,ni);\
H3.resize(ni,ni);\
H4.resize(ni,ni);\
gettimeofday(&t1, NULL);\
std::cout << s << '\t'; \
for(size_t i=0; i<NB; i++) {

#define END_BENCH \
}\
gettimeofday(&t2, NULL); \
std::cout <<  (t2.tv_sec-t1.tv_sec) + (t2.tv_usec-t1.tv_usec)/1000000. << " s" << std::endl; \
myfile << ni <<'\t' << (t2.tv_sec-t1.tv_sec) + (t2.tv_usec-t1.tv_usec)/1000000. << std::endl;\
}\
myfile.close();\

#ifdef MLCPP
#define DEFV(vi) \
cmatrix H##vi (N, N); \
for (size_t i=0; i<N; i++) \
for (size_t j=0; j<N; j++) \
H##vi(i,j) = 1;
#else
#define DEFV(vi) \
MatrixXcd H##vi (N, N); \
for (size_t i=0; i<N; i++) \
for (size_t j=0; j<N; j++) \
H##vi(i,j) = 1;
#endif


#ifdef MLCPP
using namespace mlcpp;
#else
using namespace Eigen;
#endif


#endif
