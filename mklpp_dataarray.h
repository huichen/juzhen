#ifndef MKLPP_DATAARRAY_H
#define MKLPP_DATAARRAY_H
#include <string.h>
#include <assert.h>
#include <memory>

namespace mklpp {

/////////////////////////////////////////////////////////////////////////////
/* DataArray
   wrapper for auto_ptr */
template<typename DataType> class DataArray {
public:
  DataArray(size_t s) {
    _Data = (DataType *) malloc(s*sizeof(DataType));
    assert (_Data);
    size = s;
  }

  template<typename T> DataArray(const DataArray<T> &da) {
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

  DataArray(const DataArray<DataType> &da) {
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

  template<typename T> DataArray(const DataArray<T> &da, const size_t s) {
    _Data = (DataType *) malloc(s*sizeof(DataType));
    assert (_Data);
    if (&da) {
      size = (s < da.size? s : da.size);
      for(size_t i=0; i<size; i++) _Data[i] = da[i]; 
    }
  }

  DataArray(const DataArray<DataType> &da, const size_t s) {
    _Data = (DataType *) malloc(s*sizeof(DataType));
    assert (_Data);
    if (&da) {
      size = (s < da.size? s : da.size);
      memcpy(_Data, da._Data, size*sizeof(DataType));
    }
  }

  template<typename T> DataArray(const T *da, const size_t s) {
    _Data = (DataType *) malloc(s*sizeof(DataType));
    assert (_Data);
    if (da) for (size_t i=0; i<s; i++) _Data[i] = da[i];
    size = s;
  }

  DataArray(const DataType *da, const size_t s) {
    _Data = (DataType *) malloc(s*sizeof(DataType));
    assert (_Data);
    if (da) memcpy(_Data, da, s*sizeof(DataType));
    size = s;
  }
   
  ~DataArray() {
    if (_Data) free(_Data);
  } 

  DataType& operator[](const size_t i) {
    return _Data[i];
  }
 
  DataType& operator[](const size_t i) const {
    return _Data[i];
  }

  DataType * _Data;
  size_t size;
};

template<typename DataType>
class DataPtr{
public:
  typedef std::auto_ptr<DataArray<DataType> > Type;
};
/////////////////////////////////////////////////////////////////////////////

}
#endif
