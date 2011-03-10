/*
+---------------------------------------------------------------------------+
|  Juzhen: C++ library for linear algebra                                           |
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
#ifndef BENCHMARK_BENCHMARK_H_
#define BENCHMARK_BENCHMARK_H_

#include <sys/time.h>
#include <stdio.h>

#ifdef MLCPP
#include <juzhen.h>
#else
#include <Eigen/Dense>
#endif

#include <string>

#define NB 10

const int num_digits = 2;

static int counter = 0;

/**
 * Convert an integer to a string and fill all prefix spaces with zero.
 * Example:
 * FixedLengthInteger(3, 1) returns "001"
 * FixedLengthInteger(3, 34) returns "034"
 */
std::string FixedLengthInteger(int num_digits, int input) {
  std::string target = "";
  for (int i = 0; i < num_digits; i++) {
    target = static_cast<char>(input % 10 + '0') + target;
    input = input / 10;
  }
  return target;
}

#ifdef MLCPP
#define FNAME "juzhen_"
#else
#define FNAME "eigen_"
#endif

#ifdef MLCPP
#define RESIZE_MATRIX M1.Resize(ni, ni); \
M2.Resize(ni, ni); \
M3.Resize(ni, ni); \
M4.Resize(ni, ni);
#else
#define RESIZE_MATRIX M1.resize(ni, ni); \
M2.resize(ni, ni); \
M3.resize(ni, ni); \
M4.resize(ni, ni);
#endif

#define BEGIN_BENCH(s) \
{\
{\
char out[100];\
snprintf(out, sizeof(out), "%s%s", PREFIX, \
    FixedLengthInteger(num_digits, counter).c_str());\
snprintf(filename, sizeof(filename), "%s.plt", out);\
myfile = fopen(filename, "w");\
char output_string[1000];\
snprintf(output_string, sizeof(output_string), \
"set term png\n"\
"set out \"%s.png\"\n"\
"set title \"%s\"\n"\
"set xlabel \"N\"\n"\
"set ylabel \"Seconds\"\n"\
"plot 'juzhen_%s.txt' using 1:2 title 'Juzhen' w l, "\
"'eigen_%s.txt' using 1:2 title 'Eigen' w l ", \
    out, s, out, out);\
fprintf(myfile, "%s", output_string);\
fclose(myfile);\
}\
{\
char out[100];\
snprintf(out, sizeof(out), "%s%s", PREFIX, \
    FixedLengthInteger(num_digits, counter++).c_str());\
snprintf(filename, sizeof(filename), "%s%s.txt", FNAME, out);\
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
printf("N = %ld\t", ni);\
printf("%f s\n", (t2.tv_sec-t1.tv_sec) + (t2.tv_usec-t1.tv_usec)/1000000.);\
fprintf(myfile, "%zu\t%f\n", ni, (t2.tv_sec-t1.tv_sec) \
                                 + (t2.tv_usec-t1.tv_usec)/1000000.);\
}\
fclose(myfile);\
}

#ifdef MLCPP
#define DEFV(vi) \
zmatrix M##vi(MAXN, MAXN); \
for (size_t i = 0; i < MAXN; i++) \
for (size_t j = 0; j < MAXN; j++) \
M##vi(i, j) = 1;
#else
#define DEFV(vi) \
MatrixXcd M##vi(MAXN, MAXN); \
for (size_t i = 0; i < MAXN; i++) \
for (size_t j = 0; j < MAXN; j++) \
M##vi(i, j) = 1;
#endif


#ifdef MLCPP
using juzhen::zmatrix;
#else
using Eigen::MatrixXcd;
#endif

#endif  // BENCHMARK_BENCHMARK_H_
