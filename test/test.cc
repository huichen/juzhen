/*
+---------------------------------------------------------------------------+
|  Juzhen: C++ library for linear algebra                                   |
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

#include <juzhen.h>

using juzhen::CS;
using juzhen::CD;

using juzhen::smatrix;
using juzhen::dmatrix;
using juzhen::cmatrix;
using juzhen::zmatrix;

using juzhen::identity_smatrix;
using juzhen::identity_dmatrix;
using juzhen::identity_cmatrix;
using juzhen::identity_zmatrix;

using juzhen::svector;
using juzhen::dvector;
using juzhen::cvector;
using juzhen::zvector;

int main(int argc, char *argv[]) {

#include "test_inc.h"

  RUN_TEST
  return 0;
}
