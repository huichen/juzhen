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

#ifndef MLCPP_TEST_H_  // NOLINT
#define MLCPP_TEST_H_
#include <iostream>  // NOLINT stream is a much more convenient choice here
#include <vector>
#include <string>

/**
 * test if a matrix lhs has print form rhs
*/
#define V(lhs, rhs) \
  std::cout << "Test case " << testcount++ <<":";\
  out.str(""); \
  out << lhs; \
  if (out.str() != std::string(rhs)) { \
  UnitTest::num_fail_++; \
  return_value = false;\
  std::cout << " FAILED!" << std::endl;\
  std::cout << "  Expect " << "[" << Serialize(rhs) << "]" << std::endl; \
  std::cout << "  Got    " << "[" << Serialize(String(lhs)) << "]" \
  << std::endl; \
  } else { \
  UnitTest::num_success_++; \
  std::cout << " passed." << std::endl; \
  }

/**
 * test if a matrix lhs equals to matrix rhs
*/
#define VME(lhs, rhs) \
std::cout << "Test case " << testcount++ <<":";\
if (lhs != rhs) { \
UnitTest::num_fail_++; \
return_value = false;\
std::cout << " FAILED!" << std::endl;\
  std::cout << "  Expect " << "[" << Serialize(String(lhs)) << "]" \
            << std::endl; \
  std::cout << "       = ["<< Serialize(String(rhs)) << "]"\
            << std::endl; \
} else { \
UnitTest::num_success_++; \
std::cout << " passed." << std::endl; \
}

/**
 * test if double lhs equals to double rhs
*/
#define VDE(lhs, rhs) \
std::cout << "Test case " << testcount++ <<":";\
if (lhs-rhs < -1e-4 || lhs-rhs > 1e-4) { \
UnitTest::num_fail_++; \
return_value = false;\
std::cout << " FAILED!" << std::endl;\
  std::cout << "  Expect " << lhs << " = " << rhs <<std::endl ; \
} else { \
UnitTest::num_success_++; \
std::cout << " passed." << std::endl; \
}

/**
 * test if two floats equal
*/
#define VSE(lhs, rhs) \
std::cout << "Test case " << testcount++ <<":";\
if (lhs-rhs < -1e-4 || lhs-rhs > 1e-4) { \
UnitTest::num_fail_++; \
return_value = false;\
std::cout << " FAILED!" << std::endl;\
  std::cout << "  Expect " << lhs << " = " << rhs <<std::endl ; \
} else { \
UnitTest::num_success_++; \
std::cout << " passed." << std::endl; \
}

/**
 * test if bool lhs is true
*/
#define VB(lhs) \
std::cout << "Test case " << testcount++ <<":";\
if (!(lhs)) { \
UnitTest::num_fail_++; \
return_value = false;\
std::cout << " FAILED!" << std::endl;\
  std::cout << "  Expect true." << std::endl; \
} else { \
UnitTest::num_success_++; \
std::cout << " passed." << std::endl; \
}

/**
 * all tests must begin with BEGIN_TEST
 * test is class name, classname is the class' string name
*/
#define BEGIN_TEST(test, classname)  \
class test : public TestFunction { \
  public: \
  std::string m_name; \
  test() : m_name(classname) {} \
  virtual void Run() { \
  bool return_value = true; \
  std::cout << \
    "======================================================================" \
    << std::endl; \
  std::cout << "Begin testing " << m_name << std::endl; \
  int testcount = 1; \
  std::ostringstream out;

/**
 * all tests must end with END_TEST
*/
#define END_TEST(test)  \
std::cout << m_name << (return_value?" passed.":" FAILED!") << std::endl; \
} \
}; \
static test s_TestFunction_##test; \
UnitTest::AddTest(&s_TestFunction_##test);

#define RUN_TEST UnitTest::Run();

/////////////////////////////////////////////////////////////////////////////
/**
 * Functor class. All unit tests must be derived from it.
 */
class TestFunction {
 public:
  virtual ~TestFunction() {}
  /**
   * Unit test needs to overload this virtual function, which is defined by
   * the macros BEGIN_TEST and END_TEST automatically.
   */
  virtual void Run() {}
};
/////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////
/**
 * UnitTest maintains all unit tests.
 */
class UnitTest {
 public:
  /**
   * Run all unit tests.
   */
  static void Run();

  /**
   * Add a test.
   */
  static void AddTest(TestFunction *t);

  /**
   * Number of failed tests.
   */
  static int num_fail_;

  /**
   * Number of successful tests.
   */
  static int num_success_;

 private:

  /**
   * Vector that contains all tests.
   */
  static std::vector<TestFunction*> tests_;
};

void UnitTest::Run() {
  for (size_t i = 0; i < UnitTest::tests_.size(); i++)
    UnitTest::tests_[i]->Run();
  std::cout <<
    "======================================================================"
    << std::endl; \
  std::cout << "Total unit tests: "
            << UnitTest::num_success_ + UnitTest::num_fail_
            << std::endl;
  std::cout << "Successful: " << UnitTest::num_success_ << std::endl;
  std::cout << "Failed: " << UnitTest::num_fail_ << std::endl;
}

void UnitTest::AddTest(TestFunction *t) {
  UnitTest::tests_.push_back(t);
}

/**
 * initialize static member
*/
std::vector<TestFunction*> UnitTest::tests_;
int UnitTest::num_fail_ = 0;
int UnitTest::num_success_ = 0;

/////////////////////////////////////////////////////////////////////////////

/**
 * replacing '\n' with "\\n" in a string for easy readable output
*/
std::string Serialize(std::string s) {
  std::string os;
  for (size_t i = 0; i < s.size(); i++)
    if (s[i] == '\n')
      os.append("\\n");
    else
      os.append(s.begin()+i, s.begin()+i+1);
  return os;
}
#endif  // MLCPP_TEST_H_  // NOLINT
