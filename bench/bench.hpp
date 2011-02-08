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

#define MAXN 1400
#define MINN 100
#define STEPN 100 
#define NB 10

static int counter=0;

#ifdef MLCPP 
#define FNAME "mlcpp_"
#else
#define FNAME "eigen_"
#endif

#define BEGIN_BENCH(s) \
{\
  DEFV(1)\
  DEFV(2)\
  DEFV(3)\
  DEFV(4)\
{\
std::stringstream out;\
out << counter;\
filename =out.str()+ ".plt";\
myfile.open(filename.c_str());\
myfile << "set term png \nset out \"" + out.str()+".png\"\nset title \"" + s + "\"\nset xlabel \"N\" \nset ylabel \"Seconds\" \nplot 'mlcpp_"+out.str() + ".txt' using 1:2 title 'mlcpp' w l, 'eigen_"+out.str() + ".txt' using 1:2 title 'eigen' w l ";\
myfile.close();\
}\
{\
std::stringstream out;\
out << counter++;\
filename =FNAME+out.str()+ ".txt";\
}\
myfile.open(filename.c_str());\
std::cout << "==============================" << std::endl;\
for (size_t ni=MINN; ni<=MAXN; ni+=STEPN) {\
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
}

#ifdef MLCPP
#define DEFV(vi) \
cmatrix H##vi (MAXN, MAXN); \
for (size_t i=0; i<MAXN; i++) \
for (size_t j=0; j<MAXN; j++) \
H##vi(i,j) = 1;
#else
#define DEFV(vi) \
MatrixXcd H##vi (MAXN, MAXN); \
for (size_t i=0; i<MAXN; i++) \
for (size_t j=0; j<MAXN; j++) \
H##vi(i,j) = 1;
#endif


#ifdef MLCPP
using namespace mlcpp;
#else
using namespace Eigen;
#endif


#endif
