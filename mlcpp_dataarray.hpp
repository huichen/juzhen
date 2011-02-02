#ifndef MLCPP_DATAARRAY_HPP
#define MLCPP_DATAARRAY_HPP
#include <string.h>
#include <assert.h>
#include <memory>
#include <stdlib.h>

namespace mlcpp {

/////////////////////////////////////////////////////////////////////////////
/** 
 * DataArray is a wrapper for raw array. The resources are stored in and 
 * handled by auto_ptr to prevent memory leaks.
 */
template<typename DataType> 
class DataArray {
public:
  /**
   * Default constructor.
   */
  DataArray();

  /**
   * Construct an array of size s. 
   */
  DataArray(size_t s);
 
  /**
   * Construct an array from another array.  
   */
  template<typename T> 
  DataArray(const DataArray<T> &da);
 
  /**
   * Construct an array from another array. 
   */
  DataArray(const DataArray<DataType> &da);
 
  /**
   * Construct an array from another array. The size of the new array
   * is the larger of da's size and s. 
   */
  template<typename T> 
  DataArray(const DataArray<T> &da, const size_t s);
 
  /**
   * Construct an array from another array. The size of the new array
   * is the larger of da's size and s. 
   */
  DataArray(const DataArray<DataType> &da, const size_t s);
 
  /**
   * Construct an array from a raw array. The size of the new array
   * is the larger of the raw array's size and s. 
   */
  template<typename T> 
  DataArray(const T *da, const size_t s);

  /**
   * Construct an array from a raw array. The size of the new array
   * is the larger of the raw array's size and s. 
   */
  DataArray(const DataType *da, const size_t s);
   
  /**
   * Destructor frees allocated memory.
   */
  ~DataArray();
 
  /**
   * Returns the reference of ith index in the raw array.
   */
  DataType& operator[](const size_t i); 
 
  /**
   * Returns the reference of ith index in the raw array.
   */
  DataType& operator[](const size_t i) const;

  /**
   * Returns the pointer of raw array. This is pretty useful when calling 
   * low level blas/lapack functions.
   */
  DataType * getDataPtr();

private:
  /**
   * Pointer to the raw array.
   */
  DataType * _Data;

  /**
   * Length of the raw array. It's a actual number of total items in the array.
   */
  size_t size;
};

template<typename DataType>
DataArray<DataType>::DataArray() {
  _Data = NULL;
  size = 0;
}

template<typename DataType>
DataArray<DataType>::DataArray(size_t s) {
  _Data = (DataType *) malloc(s*sizeof(DataType));
  assert (_Data);
  size = s;
}

template<typename DataType>
template<typename T> 
DataArray<DataType>::DataArray(const DataArray<T> &da) {
  if (&da && da.size >0) {
    _Data = (DataType *) malloc(da.size*sizeof(DataType));
    assert (_Data);
    size = da.size;
    if (da._Data) for(size_t i=0; i< da.size; i++) _Data[i] = da[i]; 
  } else {
    size = 0;
    _Data = NULL;
  }
}

template<typename DataType>
DataArray<DataType>::DataArray(const DataArray<DataType> &da) {
  if (&da && da.size >0) {
    _Data = (DataType *) malloc(da.size*sizeof(DataType));
    assert (_Data);
    size = da.size;
    if (da._Data) memcpy(_Data, da._Data, da.size*sizeof(DataType));
  } else {
    size = 0;
    _Data = NULL;
  }
}

template<typename DataType>
template<typename T> 
DataArray<DataType>::DataArray(const DataArray<T> &da, const size_t s) {
  if ((&da && da.size>0)||s>0) {
    size_t reals;
    if (&da)  reals = da.size>s?da.size:s;
    else reals = s;
    _Data = (DataType *) malloc(reals*sizeof(DataType));
    assert (_Data);
    size = reals;
    if (&da && da._Data) for(size_t i=0; i< da.size; i++) _Data[i] = da[i]; 
  } else {
    size = 0;
    _Data = NULL;
  }
}

template<typename DataType>
DataArray<DataType>::DataArray(const DataArray<DataType> &da, 
                               const size_t s) {
  if ((&da && da.size>0)||s>0) {
    size_t reals;
    if (&da)  reals = da.size>s?da.size:s;
    else reals = s;
    _Data = (DataType *) malloc(reals*sizeof(DataType));
    assert (_Data);
    size = reals;
    if (&da && da._Data) memcpy(_Data, da._Data, da.size*sizeof(DataType));
  } else {
    size = 0;
    _Data = NULL;
  }
}

template<typename DataType>
template<typename T> 
DataArray<DataType>::DataArray(const T *da, const size_t s) {
  if (s > 0) {
    _Data = (DataType *) malloc(s*sizeof(DataType));
    assert (_Data);
    if (da) for (size_t i=0; i<s; i++) _Data[i] = da[i];
    size = s;
  } else { 
    size = 0;
    _Data = NULL;
  }
}

template<typename DataType>
DataArray<DataType>::DataArray(const DataType *da, const size_t s) {
  if (s > 0) {
    _Data = (DataType *) malloc(s*sizeof(DataType));
    assert (_Data);
    if (da) memcpy(_Data, da, s*sizeof(DataType));
    size = s;
  } else { 
    size = 0;
    _Data = NULL;
  } 
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

template<typename DataType>
DataType * DataArray<DataType>::getDataPtr() {
  return _Data;
}
/////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////
/**
 *  DataPtr class contains the type for DataArray's auto pointer.
 */
template<typename DataType>
class DataPtr{
public:
  /**
   *  Type is the type for DataArray's auto pointer.
   */
  typedef std::auto_ptr<DataArray<DataType> > Type;
};
/////////////////////////////////////////////////////////////////////////////

}
#endif
