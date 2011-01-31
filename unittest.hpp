#ifndef UNITTEST_HPP
#define UNITTEST_HPP
#include <vector>
#include <typeinfo>
#include <string>

#define VE(lhs, rhs) \
std::cout << "Test case " << testcount++ <<":";\
if (lhs!=rhs) { \
UnitTest::m_nfail++; \
return_value = false;\
std::cout << " FAILED!" << std::endl;\
  std::cout << "  Expect " << "[" << serialize(toString(lhs)) << "]" << std::endl; \ 
  std::cout << "       = ["<< serialize(toString(rhs)) <<"]" <<std::endl ; \
} else { \
UnitTest::m_nsucc++; \
std::cout << " passed" << std::endl; \
}

#define V(lhs, rhs) \
  std::cout << "Test case " << testcount++ <<":";\
  out.str(""); \
  out << lhs; \
  if (out.str()!=std::string(rhs)) { \
  UnitTest::m_nfail++; \
  return_value = false;\
  std::cout << " FAILED!" << std::endl;\
  std::cout << "  Expect " << "[" << serialize(rhs) << "]" << std::endl; \
  std::cout << "  Got    " << "[" << serialize(toString(lhs)) << "]" << std::endl; \
  } else { \
  UnitTest::m_nsucc++; \
  std::cout << " passed." << std::endl; \
  }


#define BEGIN_TEST(test, classname)  \
class test : public TestFunction { \
public: \
  std::string m_name; \
  test () : m_name(classname) {} \
  virtual void Run() { \
  bool return_value = true; \
  std::cout << "======================================================================" << std::endl; \
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
std::vector<TestFunction*> UnitTest::m_tests;\
int UnitTest::m_nfail = 0; \
int UnitTest::m_nsucc = 0; 

/////////////////////////////////////////////////////////////////////////////
/* TestFunction class */
class TestFunction {
public:
  virtual void Run() {}
};
/////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////
/* UnitTest class */
class UnitTest {
public:
  void Run();
  static void AddTest(TestFunction *t); 
  static int m_nfail;
  static int m_nsucc;

private:
  static std::vector<TestFunction*> m_tests; 
};

void UnitTest::Run() {
  for(size_t i=0; i<UnitTest::m_tests.size(); i++)
    UnitTest::m_tests[i]->Run();
  std::cout << "======================================================================" << std::endl; \
  std::cout << "Total unit tests: " << UnitTest::m_nsucc + UnitTest::m_nfail << std::endl; 
  std::cout << "Successful: " << UnitTest::m_nsucc<< std::endl; 
  std::cout << "Failed: " << UnitTest::m_nfail << std::endl; 
} 

void UnitTest::AddTest(TestFunction *t) {
  UnitTest::m_tests.push_back(t);
}
/////////////////////////////////////////////////////////////////////////////


std::string serialize(std::string s) {
  std::string os;
  for (size_t i=0; i<s.size(); i++)
    if (s[i]=='\n') os.append("\\n");
    else os.append(s.begin()+i, s.begin()+i+1); 
  return os;
}

#endif
