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

#include "test.h"

#include <iostream>
#include <vector>
#include <sstream>

#include <mlcpp.h>

using mlcpp::CS;
using mlcpp::CD;

using mlcpp::smatrix;
using mlcpp::dmatrix;
using mlcpp::cmatrix;
using mlcpp::zmatrix;

using mlcpp::identity_smatrix;
using mlcpp::identity_dmatrix;
using mlcpp::identity_cmatrix;
using mlcpp::identity_zmatrix;

using mlcpp::svector;
using mlcpp::dvector;
using mlcpp::cvector;
using mlcpp::zvector;

int main() {

#include "test_inc.h"

  RUN_TEST
  return 0;
}
