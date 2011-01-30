#ifndef MKLPP_DATAARRAY_HPP
#define MKLPP_DATAARRAY_HPP
#include <string.h>
#include <assert.h>
#include <memory>

namespace mklpp {

/////////////////////////////////////////////////////////////////////////////
/* DataArray
   wrapper for auto_ptr */
template<typename DataType> 
class DataArray {
public:
  DataType * _Data;
  size_t size;

  DataArray(size_t s);
 
  template<typename T> 
  DataArray(const DataArray<T> &da);
 
  DataArray(const DataArray<DataType> &da);
 
  template<typename T> 
  DataArray(const DataArray<T> &da, const size_t s);
 
  DataArray(const DataArray<DataType> &da, const size_t s);
 
  template<typename T> 
  DataArray(const T *da, const size_t s);

  DataArray(const DataType *da, const size_t s);
   
  ~DataArray();
 
  DataType& operator[](const size_t i); 
 
  DataType& operator[](const size_t i) const;
};

template<typename DataType>
DataArray<DataType>::DataArray(size_t s) {
  _Data = (DataType *) malloc(s*sizeof(DataType));
  assert (_Data);
  size = s;
}

template<typename DataType>
template<typename T> 
DataArray<DataType>::DataArray(const DataArray<T> &da) {
  if (&da) {
    _Data = (DataType *) malloc(da.size*sizeof(DataType));
    assert (_Data);
    if (da._Data) for(size_t i=0; i<da.size; i++) _Data[i] = da[i]; 
    size = da.size;
  } else {
    size = 0;
    _Data = NULL;
  }
}

template<typename DataType>
DataArray<DataType>::DataArray(const DataArray<DataType> &da) {
  if (&da) {
    _Data = (DataType *) malloc(da.size*sizeof(DataType));
    assert (_Data);
    if (da._Data) memcpy(_Data, da._Data, da.size*sizeof(DataType));
    size = da.size;
  } else {
    size = 0;
    _Data = NULL;
  }
}

template<typename DataType>
template<typename T> 
DataArray<DataType>::DataArray(const DataArray<T> &da, const size_t s) {
  _Data = (DataType *) malloc(s*sizeof(DataType));
  assert (_Data);
  if (&da) {
    size = (s < da.size? s : da.size);
    for(size_t i=0; i<size; i++) _Data[i] = da[i]; 
  }
}

template<typename DataType>
DataArray<DataType>::DataArray(const DataArray<DataType> &da, 
                               const size_t s) {
  _Data = (DataType *) malloc(s*sizeof(DataType));
  assert (_Data);
  if (&da) {
    size = (s < da.size? s : da.size);
    memcpy(_Data, da._Data, size*sizeof(DataType));
  }
}

template<typename DataType>
template<typename T> 
DataArray<DataType>::DataArray(const T *da, const size_t s) {
  _Data = (DataType *) malloc(s*sizeof(DataType));
  assert (_Data);
  if (da) for (size_t i=0; i<s; i++) _Data[i] = da[i];
  size = s;
}

template<typename DataType>
DataArray<DataType>::DataArray(const DataType *da, const size_t s) {
  _Data = (DataType *) malloc(s*sizeof(DataType));
  assert (_Data);
  if (da) memcpy(_Data, da, s*sizeof(DataType));
  size = s;
}
 
template<typename DataType>
DataArray<DataType>::~DataArray() {
  if (_Data) free(_Data);
} 

template<typename DataType>
DataType& DataArray<DataType>::operator[](const size_t i) {
  return _Data[i];
}
 
template<typename DataType>
DataType& DataArray<DataType>::operator[](const size_t i) const {
  return _Data[i];
}
/////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////
/* DataPtr class*/
template<typename DataType>
class DataPtr{
public:
  typedef std::auto_ptr<DataArray<DataType> > Type;
};
/////////////////////////////////////////////////////////////////////////////

}
#endif
