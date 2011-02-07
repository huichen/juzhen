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

  /**
   * Pointer to the raw array.
   */
  DataType * m_data;

  /**
   * Length of the raw array. It's a actual number of total items in the array.
   */
  size_t m_size;
};

template<typename DataType>
DataArray<DataType>::DataArray() {
  m_data = NULL;
  m_size = 0;
}

template<typename DataType>
DataArray<DataType>::DataArray(size_t s) {
  m_data = (DataType *) malloc(s*sizeof(DataType));
#ifdef PRINT_MALLOC
  std::cout << "Malloc is called." << std::endl;
#endif
  assert (m_data);
  m_size = s;
}

template<typename DataType>
template<typename T> 
DataArray<DataType>::DataArray(const DataArray<T> &da) {
  if (&da && da.m_size >0) {
    m_data = (DataType *) malloc(da.m_size*sizeof(DataType));
#ifdef PRINT_MALLOC
  std::cout << "Malloc is called." << std::endl;
#endif
    assert (m_data);
    m_size = da.m_size;
    if (da.m_data) {
      DataType *p1 = m_data;
      T *p2 = da.m_data;
      size_t s = da.m_size;
      for(size_t i=0; i< s; i++) *(p1++) = *(p2++); 
    }
  } else {
    m_size = 0;
    m_data = NULL;
  }
}

template<typename DataType>
DataArray<DataType>::DataArray(const DataArray<DataType> &da) {
  if (&da && da.m_size >0) {
    m_data = (DataType *) malloc(da.m_size*sizeof(DataType));
#ifdef PRINT_MALLOC
  std::cout << "Malloc is called." << std::endl;
#endif
    assert (m_data);
    m_size = da.m_size;
    if (da.m_data) memcpy(m_data, da.m_data, da.m_size*sizeof(DataType));
  } else {
    m_size = 0;
    m_data = NULL;
  }
}

template<typename DataType>
template<typename T> 
DataArray<DataType>::DataArray(const DataArray<T> &da, const size_t s) {
  if ((&da && da.m_size>0)||s>0) {
    size_t reals;
    if (&da)  reals = da.m_size>s?da.m_size:s;
    else reals = s;
    m_data = (DataType *) malloc(reals*sizeof(DataType));
#ifdef PRINT_MALLOC
  std::cout << "Malloc is called." << std::endl;
#endif
    assert (m_data);
    m_size = reals;
    if (&da && da.m_data) {
      DataType *p1 = m_data;
      T *p2 = da.m_data;
      size_t s1 = da.m_size;
      for(size_t i=0; i< s1; i++) *(p1++) = *(p2++); 
    }
  } else {
    m_size = 0;
    m_data = NULL;
  }
}

template<typename DataType>
DataArray<DataType>::DataArray(const DataArray<DataType> &da, 
                               const size_t s) {
  if ((&da && da.m_size>0)||s>0) {
    size_t reals;
    if (&da)  reals = da.m_size>s?da.m_size:s;
    else reals = s;
    m_data = (DataType *) malloc(reals*sizeof(DataType));
#ifdef PRINT_MALLOC
  std::cout << "Malloc is called." << std::endl;
#endif
    assert (m_data);
    m_size = reals;
    if (&da && da.m_data) memcpy(m_data, da.m_data, da.m_size*sizeof(DataType));
  } else {
    m_size = 0;
    m_data = NULL;
  }
}

template<typename DataType>
template<typename T> 
DataArray<DataType>::DataArray(const T *da, const size_t s) {
  if (s > 0) {
    m_data = (DataType *) malloc(s*sizeof(DataType));
#ifdef PRINT_MALLOC
  std::cout << "Malloc is called." << std::endl;
#endif
    assert (m_data);
    if (da) {
      DataType *p1 = m_data;
      const T *p2 = da;
      for(size_t i=0; i<s; i++) *(p1++) = *(p2++); 
    }
    m_size = s;
  } else { 
    m_size = 0;
    m_data = NULL;
  }
}

template<typename DataType>
DataArray<DataType>::DataArray(const DataType *da, const size_t s) {
  if (s > 0) {
    m_data = (DataType *) malloc(s*sizeof(DataType));
#ifdef PRINT_MALLOC
  std::cout << "Malloc is called." << std::endl;
#endif
    assert (m_data);
    if (da) memcpy(m_data, da, s*sizeof(DataType));
    m_size = s;
  } else { 
    m_size = 0;
    m_data = NULL;
  } 
}
 
template<typename DataType>
DataArray<DataType>::~DataArray() {
  if (m_data) free(m_data);
} 

template<typename DataType>
DataType& DataArray<DataType>::operator[](const size_t i) {
  return m_data[i];
}
 
template<typename DataType>
DataType& DataArray<DataType>::operator[](const size_t i) const {
  return m_data[i];
}

template<typename DataType>
inline DataType * DataArray<DataType>::getDataPtr() {
  return m_data;
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
