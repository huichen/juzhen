/*
+---------------------------------------------------------------------------+
|  Matrix Library for C++ (mlcpp)                                           |
+---------------------------------------------------------------------------+
|                                                                           |
|  Copyright 2011 Hui Chen                                                  |
|                                                                           |
|  Licensed under the Apache License, Version 2.0 (the "License");          |
|  you may not use this file except in compliance with the License.         |
|  You may obtain a copy of the License at                                  |
|                                                                           |
|      http://www.apache.org/licenses/LICENSE-2.0                           |
|                                                                           |
|  Unless required by applicable law or agreed to in writing, software      |
|  distributed under the License is distributed on an "AS IS" BASIS,        |
|  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. |
|  See the License for the specific language governing permissions and      |
|  limitations under the License.                                           |
|                                                                           |
+---------------------------------------------------------------------------+
*/
#ifndef BENCH_BENCH_H_
#define BENCH_BENCH_H_

#include <sys/time.h>
#include <stdio.h>

#ifdef MLCPP
#include <mlcpp.h>
#else
#include <Eigen/Dense>
#endif

#define MAXN 100
#define MINN 2
#define STEPN 1
#define NB 10

static int counter = 0;

#define PREFIX "small_"

#ifdef MLCPP
#define FNAME "mlcpp_"
#else
#define FNAME "eigen_"
#endif

#ifdef MLCPP
#define RESIZE_MATRIX H1.Resize(ni, ni); \
H2.Resize(ni, ni); \
H3.Resize(ni, ni); \
H4.Resize(ni, ni);
#else
#define RESIZE_MATRIX H1.resize(ni, ni); \
H2.resize(ni, ni); \
H3.resize(ni, ni); \
H4.resize(ni, ni);
#endif

#define BEGIN_BENCH(s) \
{\
{\
char out[100];\
snprintf(out, sizeof(out), "%s%d", PREFIX, counter);\
snprintf(filename, sizeof(filename), "%s.plt", out);\
myfile = fopen(filename, "w");\
char output_string[1000];\
snprintf(output_string, sizeof(output_string), \
"set term png\n"\
"set out \"%s.png\"\n"\
"set title \"%s\"\n"\
"set xlabel \"N\"\n"\
"set ylabel \"Seconds\"\n"\
"plot 'mlcpp_%s.txt' using 1:2 title 'mlcpp' w l, "\
"'eigen_%s.txt' using 1:2 title 'eigen' w l ", \
    out, s, out, out);\
fprintf(myfile, "%s", output_string);\
fclose(myfile);\
}\
{\
char out[100];\
snprintf(out, sizeof(out), "%s%d", PREFIX, counter++);\
snprintf(filename, sizeof(filename), "%s%s.plt", FNAME, out);\
}\
myfile = fopen(filename, "w");\
printf("==============================\n");\
for (size_t ni = MINN; ni <= MAXN; ni += STEPN) {\
RESIZE_MATRIX; \
gettimeofday(&t1, NULL);\
printf("%s\t", s);\
for(size_t i = 0; i < NB; i++) {
#define END_BENCH \
}\
gettimeofday(&t2, NULL); \
printf("%f s\n", (t2.tv_sec-t1.tv_sec) + (t2.tv_usec-t1.tv_usec)/1000000.);\
fprintf(myfile, "%zu\t%f\n", ni, (t2.tv_sec-t1.tv_sec) \
                                 + (t2.tv_usec-t1.tv_usec)/1000000.);\
}\
fclose(myfile);\
}

#ifdef MLCPP
#define DEFV(vi) \
zmatrix H##vi(MAXN, MAXN); \
for (size_t i = 0; i < MAXN; i++) \
for (size_t j = 0; j < MAXN; j++) \
H##vi(i, j) = 1;
#else
#define DEFV(vi) \
MatrixXcd H##vi(MAXN, MAXN); \
for (size_t i = 0; i < MAXN; i++) \
for (size_t j = 0; j < MAXN; j++) \
H##vi(i, j) = 1;
#endif


#ifdef MLCPP
using mlcpp::zmatrix;
#else
using Eigen::MatrixXcd;
#endif

#endif  // BENCH_BENCH_H_
