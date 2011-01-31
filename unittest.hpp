#ifndef UNITTEST_HPP
#define UNITTEST_HPP
#include <vector>
#include <typeinfo>
#include <string>

#define VE(lhs, rhs) \
std::cout << "Test case " << testcount++ <<":";\
if (lhs!=rhs) { \
return_value = false;\
std::cout << name << " FAILED!";\
  std::cout << " Expect " << "[" << rhs << "] = ["<< rhs <<"]" <<std::endl ; \
} else { \
std::cout << " passed" << std::endl; \
}

#define V(lhs, rhs) \
  std::cout << "Test case " << testcount++ <<":";\
  out.str(""); \
  out << lhs; \
  if (out.str()!=std::string(rhs)) { \
  return_value = false;\
  std::cout << " FAILED!";\
  std::cout << " Expect " << "[" << rhs << "]"; \
  std::cout << " Got " << "[" << lhs << "]" << std::endl; \
  } else { \
  std::cout << " passed." << std::endl; \
  }


#define BEGIN_TEST(test, classname)  \
class test : public TestFunction { \
public: \
  std::string m_name; \
  test () : m_name(classname) {} \
  virtual void Run() { \
  bool return_value = true; \
  std::cout << "===================================" << std::endl; \
  std::cout << "Begin testing " << m_name << std::endl; \
  int testcount = 1; \
  std::ostringstream out; 
  

#define END_TEST(test)  \
std::cout << m_name << (return_value?" passed.":" FAILED!") << std::endl; \
} \
}; \
static test s_TestFunction_##test ; \
UnitTest::AddTest(&s_TestFunction_##test) ;

#define INIT_UNITTEST \
std::vector<TestFunction*> UnitTest::m_tests;

class TestFunction {
public:
  virtual void Run() {}
};

class UnitTest {
public:
  void Run();
  static void AddTest(TestFunction *t); 

private:
  static std::vector<TestFunction*> m_tests; 
};

void UnitTest::Run() {
  for(size_t i=0; i<UnitTest::m_tests.size(); i++)
    UnitTest::m_tests[i]->Run();
}

void UnitTest::AddTest(TestFunction *t) {
  UnitTest::m_tests.push_back(t);
}

#endif
