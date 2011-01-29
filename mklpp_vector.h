#ifndef MKLPP_VECTOR_H
#define MKLPP_VECTOR_H

#include "mklpp_matrix.h"

namespace mklpp {

template<typename DataType> class vector : public matrix<DataType> {
public:
  vector() : matrix<DataType>() { }
  vector(size_t i) : matrix<DataType>(i, 1) { }
  void resize(size_t n) { matrix<DataType>::resize(n,1);}
};

typedef vector<CD> cvector;
typedef vector<double> dvector;

template<typename DataType> std::ostream& operator<< (std::ostream& out, const vector<DataType> &m) {
  for (size_t i=0; i<m.nRow*m.nCol; i++) {
      out << m(i) << " ";
  }
  return out; 
}
}

#endif
