/** @file TypeHandling.h
 @author Graeme Lufkin (gwl@u.washington.edu)
 @date Created September 23, 2003
 @version 1.0
 */

#ifndef TYPEHANDLING_H
#define TYPEHANDLING_H

#include "Vector3D.h"

namespace TypeHandling {

enum DataTypeCode {
	int8 = 1,
	uint8,
	int16,
	uint16,
	int32,
	uint32,
	int64,
	uint64,
	float32,
	float64
};

inline std::ostream& operator<< (std::ostream& os, DataTypeCode code) {
	switch(code) {
		case int8:
			return os << "signed 8-bit integer";
		case uint8:
			return os << "unsigned 8-bit integer";
		case int16:
			return os << "signed 16-bit integer";
		case uint16:
			return os << "unsigned 16-bit integer";
		case int32:
			return os << "signed 32-bit integer";
		case uint32:
			return os << "unsigned 32-bit integer";
		case int64:
			return os << "signed 64-bit integer";
		case uint64:
			return os << "unsigned 64-bit integer";
		case float32:
			return os << "floating point (single precision)";
		case float64:
			return os << "floating point (double precision)";
		default:
			return os << "un-coded type!";
	}
}

/** This template structure allows you to get the enumeration value when
 you know the type.  To get the value for some type T, use Type2Code<T>::code
 T must be one of the defined types, or you will get an 
 'incomplete type' error from the compiler. */
template <typename T>
struct Type2Code;

template <>
struct Type2Code<char> {
	const static DataTypeCode code = int8;
};

template <>
struct Type2Code<unsigned char> {
	const static DataTypeCode code = uint8;
};

template <>
struct Type2Code<short> {
	const static DataTypeCode code = int16;
};

template <>
struct Type2Code<unsigned short> {
	const static DataTypeCode code = uint16;
};

template <>
struct Type2Code<int> {
	const static DataTypeCode code = int32;
};

template <>
struct Type2Code<unsigned int> {
	const static DataTypeCode code = uint32;
};

template <>
struct Type2Code<int64_t> {
	const static DataTypeCode code = int64;
};

template <>
struct Type2Code<u_int64_t> {
	const static DataTypeCode code = uint64;
};

template <>
struct Type2Code<float> {
	const static DataTypeCode code = float32;
	const static int bob = 34;
};

template <>
struct Type2Code<double> {
	const static DataTypeCode code = float64;
};

template <typename T>
struct Type2Code<Vector3D<T> > {
	const static DataTypeCode code = Type2Code<T>::code;
};

/// Dimensionality of a type.  Currently only scalar and vector are supported.
template <typename T>
struct Type2Dimensions {
	const static unsigned int dimensions = 1;
};

template <typename T>
struct Type2Dimensions<Vector3D<T> > {
	const static unsigned int dimensions = 3;
};

/// A TypeInformation object holds the dimensionality and type together
class TypeInformation {
public:
	unsigned int dimensions;
	DataTypeCode code;
	
	TypeInformation() : dimensions(0) { }
	
	template <typename T>
	void initialize() {
		dimensions = Type2Dimensions<T>::dimensions;
		code = Type2Code<T>::code;
	}
};

template <typename T>
inline bool checkType(const TypeInformation& info) {
	return (Type2Dimensions<T>::dimensions == info.dimensions && Type2Code<T>::code == info.code);
}

template <typename T>
inline void* allocateArray(const unsigned int dimensions, const u_int64_t N) {
	if(dimensions == 1)
		return new T[N];
	else if(dimensions == 3)
		return new Vector3D<T>[N];
	else
		//throw UnsupportedDimensionalityException;
		return 0;
}

inline void* allocateArray(const TypeInformation& info, const u_int64_t N) {
	switch(info.code) {
		case int8:
			return allocateArray<char>(info.dimensions, N);
		case uint8:
			return allocateArray<unsigned char>(info.dimensions, N);
		case int16:
			return allocateArray<short>(info.dimensions, N);
		case uint16:
			return allocateArray<unsigned short>(info.dimensions, N);
		case int32:
			return allocateArray<int>(info.dimensions, N);
		case uint32:
			return allocateArray<unsigned int>(info.dimensions, N);
		case int64:
			return allocateArray<int64_t>(info.dimensions, N);
		case uint64:
			return allocateArray<u_int64_t>(info.dimensions, N);
		case float32:
			return allocateArray<float>(info.dimensions, N);
		case float64:
			return allocateArray<double>(info.dimensions, N);
		default:
			//throw UnsupportedTypeException;
			return 0;
	}
}

template <typename T>
inline void* allocateValue(const unsigned int dimensions) {
	if(dimensions == 1)
		return new T;
	else if(dimensions == 3)
		return new Vector3D<T>;
	else
		//throw UnsupportedDimensionalityException;
		return 0;
}

inline void* allocateValue(const TypeInformation& info) {
	switch(info.code) {
		case int8:
			return allocateValue<char>(info.dimensions);
		case uint8:
			return allocateValue<unsigned char>(info.dimensions);
		case int16:
			return allocateValue<short>(info.dimensions);
		case uint16:
			return allocateValue<unsigned short>(info.dimensions);
		case int32:
			return allocateValue<int>(info.dimensions);
		case uint32:
			return allocateValue<unsigned int>(info.dimensions);
		case int64:
			return allocateValue<int64_t>(info.dimensions);
		case uint64:
			return allocateValue<u_int64_t>(info.dimensions);
		case float32:
			return allocateValue<float>(info.dimensions);
		case float64:
			return allocateValue<double>(info.dimensions);
		default:
			//throw UnsupportedTypeException;
			return 0;
	}
}

template <typename T>
inline bool deleteArray(const unsigned int dimensions, void*& data) {
	if(dimensions == 1)
		delete[] reinterpret_cast<T *>(data);
	else if(dimensions == 3)
		delete[] reinterpret_cast<Vector3D<T> *>(data);
	else
		//throw UnsupportedDimensionalityException;
		return false;
	
	data = 0;
	return true;
}

inline bool deleteArray(const TypeInformation& info, void*& data) {
	switch(info.code) {
		case int8:
			return deleteArray<char>(info.dimensions, data);
		case uint8:
			return deleteArray<unsigned char>(info.dimensions, data);
		case int16:
			return deleteArray<short>(info.dimensions, data);
		case uint16:
			return deleteArray<unsigned short>(info.dimensions, data);
		case int32:
			return deleteArray<int>(info.dimensions, data);
		case uint32:
			return deleteArray<unsigned int>(info.dimensions, data);
		case int64:
			return deleteArray<int64_t>(info.dimensions, data);
		case uint64:
			return deleteArray<u_int64_t>(info.dimensions, data);
		case float32:
			return deleteArray<float>(info.dimensions, data);
		case float64:
			return deleteArray<double>(info.dimensions, data);
		default:
			//throw UnsupportedTypeException;
			return false;
	}
}

template <typename T>
inline bool deleteValue(const unsigned int dimensions, void*& data) {
	if(dimensions == 1)
		delete reinterpret_cast<T *>(data);
	else if(dimensions == 3)
		delete reinterpret_cast<Vector3D<T> *>(data);
	else
		//throw UnsupportedDimensionalityException;
		return false;
	
	data = 0;
	return true;
}

inline bool deleteValue(const TypeInformation& info, void*& data) {
	switch(info.code) {
		case int8:
			return deleteValue<char>(info.dimensions, data);
		case uint8:
			return deleteValue<unsigned char>(info.dimensions, data);
		case int16:
			return deleteValue<short>(info.dimensions, data);
		case uint16:
			return deleteValue<unsigned short>(info.dimensions, data);
		case int32:
			return deleteValue<int>(info.dimensions, data);
		case uint32:
			return deleteValue<unsigned int>(info.dimensions, data);
		case int64:
			return deleteValue<int64_t>(info.dimensions, data);
		case uint64:
			return deleteValue<u_int64_t>(info.dimensions, data);
		case float32:
			return deleteValue<float>(info.dimensions, data);
		case float64:
			return deleteValue<double>(info.dimensions, data);
		default:
			//throw UnsupportedTypeException;
			return false;
	}
}

} //close namespace TypeHandling

#endif //TYPEHANDLING_H
