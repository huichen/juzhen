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

#ifndef SRC_CORE_DATA_ARRAY_H_
#define SRC_CORE_DATA_ARRAY_H_
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <memory>

namespace juzhen {

/////////////////////////////////////////////////////////////////////////////
/**
 * DataArray is a wrapper for raw array. The resources are stored in and
 * handled by auto_ptr to prevent memory leaks.
 */
template<typename DataType>
struct DataArray {
  /**
   * Default constructor.
   */
  DataArray();

  /**
   * Construct an array of size s.
   */
  DataArray(size_t s);  // NOLINT

  /**
   * Construct an array from another array.
   */
  template<typename T>
  DataArray(const DataArray<T> &da);  // NOLINT

  /**
   * Construct an array from another array.
   */
  DataArray(const DataArray<DataType> &da);  // NOLINT

  /**
   * Construct an array from another array. The size of the new array
   * is the larger of da's size and s.
   */
  template<typename T>
  DataArray(const DataArray<T> &da, size_t s);

  /**
   * Construct an array from another array. The size of the new array
   * is the larger of da's size and s.
   */
  DataArray(const DataArray<DataType> &da, size_t s);

  /**
   * Construct an array from a raw array. The size of the new array
   * is the larger of the raw array's size and s.
   */
  template<typename T>
  DataArray(const T *da, size_t s);

  /**
   * Construct an array from a raw array. The size of the new array
   * is the larger of the raw array's size and s.
   */
  DataArray(const DataType *da, size_t s);

  /**
   * Destructor frees allocated memory.
   */
  ~DataArray();

  /**
   * Returns the reference of ith index in the raw array.
   */
  DataType &operator[](size_t i);

  /**
   * Returns the reference of ith index in the raw array.
   */
  DataType &operator[](size_t i) const;

  /**
   * Pointer to the raw array.
   */
  DataType *data_ptr;

  /**
   * Length of the raw array. It's a actual number of total items in the array.
   */
  size_t size;
};

template<typename DataType>
DataArray<DataType>::DataArray() {
  data_ptr = NULL;
  size = 0;
}

template<typename DataType>
DataArray<DataType>::DataArray(size_t s) {
  data_ptr = reinterpret_cast<DataType *>(malloc(s*sizeof(DataType)));
#ifdef PRINT_MALLOC
  std::cout << "DataArray(size_t) is called." << std::endl;
#endif
  assert(data_ptr);
  size = s;
}

template<typename DataType>
template<typename T>
DataArray<DataType>::DataArray(const DataArray<T> &da) {
  if (&da && da.size > 0) {
    data_ptr = reinterpret_cast<DataType *>(malloc(da.size*sizeof(DataType)));
#ifdef PRINT_MALLOC
  std::cout << "DataArray(const DataArray<T> &) is called." << std::endl;
#endif
    assert(data_ptr);
    size = da.size;
    if (da.data_ptr) {
      DataType *p1 = data_ptr;
      T *p2 = da.data_ptr;
      size_t s = da.size;
      for (size_t i = 0; i < s; i++) *(p1++) = *(p2++);
    }
  } else {
    size = 0;
    data_ptr = NULL;
  }
}

template<typename DataType>
DataArray<DataType>::DataArray(const DataArray<DataType> &da) {
  if (&da && da.size > 0) {
    data_ptr = reinterpret_cast<DataType *>(malloc(da.size*sizeof(DataType)));
#ifdef PRINT_MALLOC
  std::cout << "DataArray(const DataArray<DataType> &) is called." << std::endl;
#endif
    assert(data_ptr);
    size = da.size;
    if (da.data_ptr) memcpy(data_ptr, da.data_ptr, da.size*sizeof(DataType));
  } else {
    size = 0;
    data_ptr = NULL;
  }
}

template<typename DataType>
template<typename T>
DataArray<DataType>::DataArray(const DataArray<T> &da, size_t s) {
  if ((&da && da.size > 0) || s > 0) {
    size_t reals;
    if (&da)
      reals = da.size > s ? da.size : s;
    else
      reals = s;
    data_ptr = reinterpret_cast<DataType *>(malloc(reals*sizeof(DataType)));
#ifdef PRINT_MALLOC
  std::cout << "DataArray(const DataArray<T> &, size_t) is called."
            << std::endl;
#endif
    assert(data_ptr);
    size = reals;
    if (&da && da.data_ptr) {
      DataType *p1 = data_ptr;
      T *p2 = da.data_ptr;
      size_t s1 = da.size;
      for (size_t i = 0; i < s1; i++) *(p1++) = *(p2++);
    }
  } else {
    size = 0;
    data_ptr = NULL;
  }
}

template<typename DataType>
DataArray<DataType>::DataArray(const DataArray<DataType> &da,
                               size_t s) {
  if ((&da && da.size > 0) || s > 0) {
    size_t reals;
    if (&da)
      reals = da.size > s ? da.size : s;
    else
      reals = s;
    data_ptr = reinterpret_cast<DataType *>(malloc(reals*sizeof(DataType)));
#ifdef PRINT_MALLOC
  std::cout << "DataArray(const DataArray<DataType> &, size_t) is called."
            << std::endl;
#endif
    assert(data_ptr);
    size = reals;
    if (&da && da.data_ptr) memcpy(data_ptr,
                                   da.data_ptr,
                                   da.size*sizeof(DataType));
  } else {
    size = 0;
    data_ptr = NULL;
  }
}

template<typename DataType>
template<typename T>
DataArray<DataType>::DataArray(const T *da, size_t s) {
  if (s > 0) {
    data_ptr = reinterpret_cast<DataType *>(malloc(s*sizeof(DataType)));
#ifdef PRINT_MALLOC
  std::cout << "DataArray(const T *, size_t) is called." << std::endl;
#endif
    assert(data_ptr);
    if (da) {
      DataType *p1 = data_ptr;
      const T *p2 = da;
      for (size_t i = 0; i < s; i++) *(p1++) = *(p2++);
    }
    size = s;
  } else {
    size = 0;
    data_ptr = NULL;
  }
}

template<typename DataType>
DataArray<DataType>::DataArray(const DataType *da, size_t s) {
  if (s > 0) {
    data_ptr = reinterpret_cast<DataType *>(malloc(s*sizeof(DataType)));
#ifdef PRINT_MALLOC
  std::cout << "DataArray(const DataType *, size_t) is called." << std::endl;
#endif
    assert(data_ptr);
    if (da) memcpy(data_ptr, da, s*sizeof(DataType));
    size = s;
  } else {
    size = 0;
    data_ptr = NULL;
  }
}

template<typename DataType>
DataArray<DataType>::~DataArray() {
  if (data_ptr) free(data_ptr);
}

template<typename DataType>
DataType &DataArray<DataType>::operator[](size_t i) {
  return data_ptr[i];
}

template<typename DataType>
DataType &DataArray<DataType>::operator[](size_t i) const {
  return data_ptr[i];
}
/////////////////////////////////////////////////////////////////////////////
}
#endif  // SRC_CORE_DATA_ARRAY_H_
