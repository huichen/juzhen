#include <iostream>
#include <mkl_types.h>
#include <mkl_cblas.h>
#include <mkl.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

namespace MAT {
template<typename DataType> class Matrix {
public:
  size_t nCol;
  size_t nRow;

  Matrix() : nCol(0), nRow(0) {
    _Data = NULL;  
    _DataSize = 0;
  } 

  Matrix(size_t nr, size_t nc) : nCol(nc), nRow(nr) {
    _Data = NULL;  
    _DataSize = 0;
  } 

  Matrix(const DataType *data, size_t nr, size_t nc) : nCol(nc), nRow(nr) {
    nCol = nc;
    nRow = nr;
    _Data = (DataType *) malloc(nCol*nRow*sizeof(DataType));
    assert(_Data);
    memcpy(_Data, data, nCol*nRow*sizeof(DataType));    
    _DataSize = nCol * nRow;
  } 

  DataType * getDataPtr() const {
    return _Data;
  }

  size_t getDataSize() {
    return _DataSize;
  }

  Matrix(const Matrix<DataType> &m) {
    nCol = m.nCol;
    nRow = m.nRow;
    _Data = (DataType *) malloc(nCol*nRow*sizeof(DataType));
    memcpy(_Data, m.getDataPtr(), nCol*nRow*sizeof(DataType));    
    _DataSize = nCol*nRow;
  } 

  const Matrix<DataType>& operator= (const Matrix<DataType> &rhs) {
    if (&rhs==this) return *this;
    nCol = rhs.nCol;
    nRow = rhs.nRow;
    if (nCol*nRow > _DataSize) {
      free(_Data);
      _Data = (DataType *) malloc(nCol*nRow*sizeof(DataType));
      assert(_Data);
      _DataSize = nCol*nRow;
    }
    memcpy(_Data, rhs.getDataPtr(), nCol*nRow*sizeof(DataType));    
    return *this;
  } 

  const Matrix<DataType>& operator= (const DataType *rhs) {
    if (sizeof(rhs)>_DataSize) {
      if (_Data) free(_Data);
      _Data = (DataType *) malloc( sizeof(rhs)*sizeof(DataType));
      assert(_Data);
      _DataSize = sizeof(rhs); 
    }
    memcpy(_Data, rhs, sizeof(rhs)*sizeof(DataType));
    return *this;
  } 

  void setDim(size_t nr, size_t nc) {
    if (_DataSize < nc * nr) {
      DataType * newData = (DataType *) malloc( nc*nr*sizeof(DataType));
      assert(newData);
      memcpy(newData, _Data, nCol*nRow*sizeof(DataType));
      if (_Data) free(_Data);
      _Data = newData;
      _DataSize = nc*nr;
    }
    nCol = nc;
    nRow = nr;
  } 

  DataType& operator()(size_t i, size_t j) {
    assert(i<nRow && j<nCol);
    return _Data[j*nRow + i];
  }

  DataType operator()(size_t i, size_t j) const {
    assert(i<nRow && j<nCol);
    return _Data[j*nRow + i];
  }

  const Matrix<DataType> operator*(const Matrix<DataType>& rhs) {
    assert (nCol == rhs.nRow);
    Matrix<DataType> m;

    m.setDim(nRow, rhs.nCol);
    DataType temp;
    for (size_t i=0; i<nRow; i++)
      for (size_t j=0; j<rhs.nCol; j++) {
        temp = 0;
        for (size_t k=0; k<nCol; k++) 
          temp += (*this)(i,k)*rhs(k,j);
        m(i,j) = temp;
      }

    return m;
  } 
  
private:
  DataType * _Data; 
  size_t _DataSize;
};

template<typename DataType> std::ostream& operator<< (std::ostream& out, const Matrix<DataType> &m) {
  for (size_t i=0; i<m.nRow; i++) {
    for (size_t j=0; j<m.nCol; j++)
      out << m(i,j) << " ";
    if (i!=m.nRow-1) out << std::endl;
  }
  return out; 
}
 
}

using namespace MAT;
using namespace std;

int main() {

  Matrix<double> m1;

  double a[] = {1., 2., 3., 4.};
  double b[] = {1., 2., 3., 4., 5., 6.};

  m1 = a;
  m1.setDim(2, 2);
  std::cout << m1(0,1) << std::endl;

  Matrix<double> m2 = m1;
  m2(0,1) = 99.;

  Matrix<double> m3(b, 2, 3);

  cout << m1 << endl;
  cout << m2 << endl;
  cout << m1*m2 << endl;

  cout << m3 << endl;
  cout << m1*m3 << endl;

  return 0;
}
