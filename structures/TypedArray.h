/** @file TypedArray.h
 @author Graeme Lufkin (gwl@u.washington.edu)
 @date Created September 23, 2003
 @version 1.0
 */

#ifndef TYPEDARRAY_H
#define TYPEDARRAY_H

#include "TypeHandling.h"

namespace TypeHandling {

template <typename T>
struct Type2Type {
	typedef T OriginalType;
};

/** A TypedArray is the recommended way to store an array of values
 whose dimensionality and type are determined at runtime.  It provides
 methods for accessing the array of values in safe and unsafe ways,
 and offers safe ways to get the minimum and maximum values, as taken
 from the attribute file. */
class TypedArray : public TypeInformation {
	template <typename T>
	void calculateMinMax(Type2Type<T>) {
		T* minVal = reinterpret_cast<T *>(minValue);
		T* maxVal = reinterpret_cast<T *>(maxValue);
		T* array = reinterpret_cast<T *>(data);
		*minVal = array[0];
		*maxVal = array[0];
		for(u_int64_t i = 1; i < length; ++i) {
			if(array[i] < *minVal)
				*minVal = array[i];
			if(*maxVal < array[i])
				*maxVal = array[i];
		}
	}
	
	template <typename T>
	void calculateMinMaxVector(Type2Type<T>) {
		Vector3D<T>* array = reinterpret_cast<Vector3D<T> *>(data);
		OrientedBox<T> box(array[0], array[0]);
		for(u_int64_t i = 1; i < length; ++i)
			box.grow(array[i]);
		*reinterpret_cast<Vector3D<T> *>(minValue) = box.lesser_corner;
		*reinterpret_cast<Vector3D<T> *>(maxValue) = box.greater_corner;
	}
	
	template <typename T>
	void calculateMinMax_dimensions(Type2Type<T>) {
		if(dimensions == 1)
			calculateMinMax(Type2Type<T>());
		else if(dimensions == 3)
			calculateMinMaxVector(Type2Type<T>());
		//else
			//throw UnsupportedDimensionalityException;
	}
	
public:
	
	void* minValue;
	void* maxValue;
	
	u_int64_t length;
	void* data;

	TypedArray() : minValue(0), maxValue(0), length(0), data(0) { }	

	template <typename T>
	T* getArray(Type2Type<T>) {
		//if(Type2Dimensions<T>::dimensions != dimensions || Type2Code<T>::code != code)
			//throw TypeMismatchException;
		//	return 0;
		return reinterpret_cast<T *>(data);
	}

	template <typename T>
	T getMinValue(Type2Type<T>) const {
		//if(Type2Dimensions<T>::dimensions != dimensions || Type2Code<T>::code != code)
		//	throw TypeMismatchException;
		return *reinterpret_cast<const T *>(minValue);
	}

	template <typename T>
	T getMaxValue(Type2Type<T>) const {
		//if(Type2Dimensions<T>::dimensions != dimensions || Type2Code<T>::code != code)
		//	throw TypeMismatchException;
		return *reinterpret_cast<const T *>(maxValue);
	}
	
	void release() {
		deleteArray(*this, data);
		deleteValue(*this, minValue);
		deleteValue(*this, maxValue);
	}
	
	void calculateMinMax() {
		if(length != 0) {
			if(minValue == 0)
				minValue = allocateValue(*this);
			if(maxValue == 0)
				maxValue = allocateValue(*this);
			switch(code) {
				case int8:
					return calculateMinMax_dimensions(Type2Type<char>());
				case uint8:
					return calculateMinMax_dimensions(Type2Type<unsigned char>());
				case int16:
					return calculateMinMax_dimensions(Type2Type<short>());
				case uint16:
					return calculateMinMax_dimensions(Type2Type<unsigned short>());
				case int32:
					return calculateMinMax_dimensions(Type2Type<int>());
				case uint32:
					return calculateMinMax_dimensions(Type2Type<unsigned int>());
				case int64:
					return calculateMinMax_dimensions(Type2Type<int64_t>());
				case uint64:
					return calculateMinMax_dimensions(Type2Type<u_int64_t>());
				case float32:
					return calculateMinMax_dimensions(Type2Type<float>());
				case float64:
					return calculateMinMax_dimensions(Type2Type<double>());
				//default:
					//throw UnsupportedTypeException;
			}
		}
	}
};

template <typename T>
inline double getScalarValue(const unsigned int dimensions, const void* value) {
	if(dimensions == 1)
		return static_cast<double>(*reinterpret_cast<const T *>(value));
	else if(dimensions == 3)
		return static_cast<double>(reinterpret_cast<const Vector3D<T> *>(value)->length());
	else
		//throw UnsupportedDimensionalityException;
		return 0;
}

inline double getScalarMin(const TypedArray& arr) {
	switch(arr.code) {
		case int8:
			return getScalarValue<char>(arr.dimensions, arr.minValue);
		case uint8:
			return getScalarValue<unsigned char>(arr.dimensions, arr.minValue);
		case int16:
			return getScalarValue<short>(arr.dimensions, arr.minValue);
		case uint16:
			return getScalarValue<unsigned short>(arr.dimensions, arr.minValue);
		case int32:
			return getScalarValue<int>(arr.dimensions, arr.minValue);
		case uint32:
			return getScalarValue<unsigned int>(arr.dimensions, arr.minValue);
		case int64:
			return getScalarValue<int64_t>(arr.dimensions, arr.minValue);
		case uint64:
			return getScalarValue<u_int64_t>(arr.dimensions, arr.minValue);
		case float32:
			return getScalarValue<float>(arr.dimensions, arr.minValue);
		case float64:
			return getScalarValue<double>(arr.dimensions, arr.minValue);
		default:
			//throw UnsupportedTypeException;
			return HUGE_VAL;
	}
}

inline double getScalarMax(const TypedArray& arr) {
	switch(arr.code) {
		case int8:
			return getScalarValue<char>(arr.dimensions, arr.maxValue);
		case uint8:
			return getScalarValue<unsigned char>(arr.dimensions, arr.maxValue);
		case int16:
			return getScalarValue<short>(arr.dimensions, arr.maxValue);
		case uint16:
			return getScalarValue<unsigned short>(arr.dimensions, arr.maxValue);
		case int32:
			return getScalarValue<int>(arr.dimensions, arr.maxValue);
		case uint32:
			return getScalarValue<unsigned int>(arr.dimensions, arr.maxValue);
		case int64:
			return getScalarValue<int64_t>(arr.dimensions, arr.maxValue);
		case uint64:
			return getScalarValue<u_int64_t>(arr.dimensions, arr.maxValue);
		case float32:
			return getScalarValue<float>(arr.dimensions, arr.maxValue);
		case float64:
			return getScalarValue<double>(arr.dimensions, arr.maxValue);
		default:
			//throw UnsupportedTypeException;
			return -HUGE_VAL;
	}
}

} //close namespace TypeHandling

#endif //TYPEDARRAY_H
