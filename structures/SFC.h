/** @file SFC.h
 Structures, defines, functions relating to space-filling curves,
 as used to build trees for particle data.
 */

#ifndef SFC_H
#define SFC_H

#include <iostream>

#include "Vector3D.h"
#include "OrientedBox.h"

namespace SFC {

typedef u_int64_t Key;

inline void printFloatBits(float f, std::ostream& os) {
	int i = *reinterpret_cast<int *>(&f);
	if(i & (1 << 31))
		os << "1 ";
	else
		os << "0 ";
	for(int j = 30; j >= 23; --j) {
		if(i & (1 << j))
			os << "1";
		else
			os << "0";
	}
	os << " ";
	for(int j = 22; j >= 0; --j) {
		if(i & (1 << j))
			os << "1";
		else
			os << "0";
	}
}

inline void printKeyBits(Key k, std::ostream& os) {
	for(int i = 63; i >= 0; --i) {
		if(k & (static_cast<Key>(1) << i))
			os << "1";
		else
			os << "0";
	}
}

inline void printIntBits(int k, std::ostream& os) {
	for(int i = 31; i >= 0; --i) {
		if(k & (1 << i))
			os << "1";
		else
			os << "0";
	}
}

/** The very first possible key a particle can take on. */
const Key firstPossibleKey = static_cast<Key>(0);
/** The very last possible key a particle can take on. */
const Key lastPossibleKey = ~(static_cast<Key>(1) << 63);

float exchangeKey[3];

/** Given the floating point numbers for the location, construct the key. 
 The key uses 21 of the 23 bits for the floats of the x, y, and z coordinates
 of the position vector.  This process will only make sense if the position
 coordinates are in the range [1,2).  The mantissa bits are taken, and interleaved
 in xyz order to form the key.  This makes the key a position on the z-ordering
 space-filling curve. 
 The second parameter is the bounding Box used to normalize the position vectors
 between 1 and 2.
 */
inline Key makeKey() {
  //unsigned int ix = *reinterpret_cast<unsigned int *>(&exchangeKey[0]);
  //unsigned int iy = *reinterpret_cast<unsigned int *>(&exchangeKey[1]);
  //unsigned int iz = *reinterpret_cast<unsigned int *>(&exchangeKey[2]);
  unsigned int ix = (unsigned int)(exchangeKey[0]*(1<<21) - exchangeKey[0]);
  unsigned int iy = (unsigned int)(exchangeKey[1]*(1<<21) - exchangeKey[1]);
  unsigned int iz = (unsigned int)(exchangeKey[2]*(1<<21) - exchangeKey[2]);
	Key key = 0;
	for(unsigned int mask = (1 << 20); mask > 0; mask >>= 1) {
		key <<= 3;
		if(ix & mask)
			key += 4;
		if(iy & mask)
			key += 2;
		if(iz & mask)
			key += 1;
	}
	return key;
}

//Key generateKey(const Vector3D<float>& v, const OrientedBox<float>& boundingBox);

inline Key generateKey(const Vector3D<float>& v, const OrientedBox<float>& boundingBox){
  Vector3D<float> d = (v - boundingBox.lesser_corner) / (boundingBox.greater_corner - boundingBox.lesser_corner); //+ Vector3D<float>(1, 1, 1);
  exchangeKey[0] = d.x;
  exchangeKey[1] = d.y;
  exchangeKey[2] = d.z;
  return makeKey();
}

/** Given a key, create a vector of floats representing a position.
 This is almost the inverse of the makeKey() function.  Since the key
 does not use all the bits, each float generated here will have its last
 two mantissa bits set zero, regardless of the values when the key was
 originally generated.
 */
inline Vector3D<float> makeVector(Key k){
	//Vector3D<float> v(1.0, 1.0, 1.0);
	//int* ix = reinterpret_cast<int *>(&v.x);
	//int* iy = reinterpret_cast<int *>(&v.y);
	//int* iz = reinterpret_cast<int *>(&v.z);
  unsigned int ix=0, iy=0, iz=0;
	for(int mask = (1 << 0); mask <= (1 << 20); mask <<= 1) {
		if(k & 4)
			ix |= mask;
		if(k & 2)
			iy |= mask;
		if(k & 1)
			iz |= mask;
		k >>= 3;	
	}
  float x = ((float)ix) / ((1<<21)-1);
  float y = ((float)iy) / ((1<<21)-1);
  float z = ((float)iz) / ((1<<21)-1);
	Vector3D<float> v(x, y, z);
	return v;
}

template <typename T>
inline OrientedBox<T> cutBoxLeft(const OrientedBox<T>& box, const int axis) {
	OrientedBox<T> newBox = box;
	switch(axis) {
		case 0: //cut x
			newBox.greater_corner.x = (box.greater_corner.x + box.lesser_corner.x) / 2.0;
			break;
		case 1: //cut y
			newBox.greater_corner.y = (box.greater_corner.y + box.lesser_corner.y) / 2.0;				
			break;
		case 2: //cut z
			newBox.greater_corner.z = (box.greater_corner.z + box.lesser_corner.z) / 2.0;
			break;
	}
	return newBox;
}

template <typename T>
inline OrientedBox<T> cutBoxRight(const OrientedBox<T>& box, const int axis) {
	OrientedBox<T> newBox = box;
	switch(axis) {
		case 0: //cut x
			newBox.lesser_corner.x = (box.greater_corner.x + box.lesser_corner.x) / 2.0;
			break;
		case 1: //cut y
			newBox.lesser_corner.y = (box.greater_corner.y + box.lesser_corner.y) / 2.0;				
			break;
		case 2: //cut z
			newBox.lesser_corner.z = (box.greater_corner.z + box.lesser_corner.z) / 2.0;
			break;
	}
	return newBox;
}

/** Increment the floating point number until the last two bits
 of the mantissa are zero.
 */
inline float bumpLastTwoBits(float f, float direction) {
	int howmany = *reinterpret_cast<int *>(&f) & 3;
	if(direction < 0) {
		switch(howmany) {
			case 1:
				howmany = 3;
				break;
			case 3:
				howmany = 1;
				break;
		}
	}
	switch(howmany) {
		case 0:
			f = nextafterf(f, direction);
		case 1:
			f = nextafterf(f, direction);
		case 2:
			f = nextafterf(f, direction);
		case 3:
			f = nextafterf(f, direction);
	}
	return f;
}

inline void bumpBox(OrientedBox<float>& b, float direction) {
	b.greater_corner.x = bumpLastTwoBits(b.greater_corner.x, direction);
	b.greater_corner.y = bumpLastTwoBits(b.greater_corner.y, direction);
	b.greater_corner.z = bumpLastTwoBits(b.greater_corner.z, direction);
	b.lesser_corner.x = bumpLastTwoBits(b.lesser_corner.x, -direction);
	b.lesser_corner.y = bumpLastTwoBits(b.lesser_corner.y, -direction);
	b.lesser_corner.z = bumpLastTwoBits(b.lesser_corner.z, -direction);	
}

template <typename T>
inline void cubize(OrientedBox<T>& b) {
	T max = b.greater_corner.x - b.lesser_corner.x;
	if((b.greater_corner.y - b.lesser_corner.y) > max)
		max = b.greater_corner.y - b.lesser_corner.y;
	if((b.greater_corner.z - b.lesser_corner.z) > max)
		max = b.greater_corner.z - b.lesser_corner.z;
	T middle = (b.greater_corner.x + b.lesser_corner.x) / 2.0;
	b.greater_corner.x = middle + max / 2.0;
	b.lesser_corner.x = middle - max / 2.0;
	middle = (b.greater_corner.y + b.lesser_corner.y) / 2.0;
	b.greater_corner.y = middle + max / 2.0;
	b.lesser_corner.y = middle - max / 2.0;
	middle = (b.greater_corner.z + b.lesser_corner.z) / 2.0;
	b.greater_corner.z = middle + max / 2.0;
	b.lesser_corner.z = middle - max / 2.0;
	
}

} //close namespace SFC

#endif //SFC_H
