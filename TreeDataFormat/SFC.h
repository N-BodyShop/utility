/** @file SFC.h
 Structures, defines, functions relating to space-filling curves,
 as used to build trees for particle data.
 */

#ifndef SFC_H
#define SFC_H

#include <iostream>

namespace SFC {

typedef unsigned long long Key;

void printFloatBits(float f, std::ostream& os) {
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

void printKeyBits(Key k, std::ostream& os) {
	for(int i = 63; i >= 0; --i) {
		if(k & (static_cast<Key>(1) << i))
			os << "1";
		else
			os << "0";
	}
}

void printIntBits(int k, std::ostream& os) {
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

/// Out of three floats, make a Morton order (z-order) space-filling curve key
// Cannot make this function inline, or g++ with -O2 or higher generates incorrect code
Key makeKey(Vector3D<float> v) {
	unsigned int ix = *reinterpret_cast<unsigned int *>(&v.x);
	unsigned int iy = *reinterpret_cast<unsigned int *>(&v.y);
	unsigned int iz = *reinterpret_cast<unsigned int *>(&v.z);
	Key key = 0;
	for(unsigned int mask = (1 << 22); mask > 2; mask >>= 1) {
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

template <typename T, typename T2>
inline Key generateKey(const Vector3D<T>& v, const OrientedBox<T2>& boundingBox) {
	return makeKey((v - boundingBox.lesser_corner) / (boundingBox.greater_corner - boundingBox.lesser_corner) + Vector3D<float>(1, 1, 1));
}

template <typename T>
OrientedBox<T> cutBoxLeft(const OrientedBox<T>& box, const int axis) {
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
OrientedBox<T> cutBoxRight(const OrientedBox<T>& box, const int axis) {
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
float bumpLastTwoBits(float f, float direction) {
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

void bumpBox(OrientedBox<float>& b, float direction) {
	b.greater_corner.x = bumpLastTwoBits(b.greater_corner.x, direction);
	b.greater_corner.y = bumpLastTwoBits(b.greater_corner.y, direction);
	b.greater_corner.z = bumpLastTwoBits(b.greater_corner.z, direction);
	b.lesser_corner.x = bumpLastTwoBits(b.lesser_corner.x, -direction);
	b.lesser_corner.y = bumpLastTwoBits(b.lesser_corner.y, -direction);
	b.lesser_corner.z = bumpLastTwoBits(b.lesser_corner.z, -direction);	
}

void cubize(OrientedBox<float>& b) {
	float max = b.greater_corner.x - b.lesser_corner.x;
	if((b.greater_corner.y - b.lesser_corner.y) > max)
		max = b.greater_corner.y - b.lesser_corner.y;
	if((b.greater_corner.z - b.lesser_corner.z) > max)
		max = b.greater_corner.z - b.lesser_corner.z;
	float middle = (b.greater_corner.x + b.lesser_corner.x) / 2.0;
	b.greater_corner.x = middle + max / 2.0;
	b.lesser_corner.x = middle - max / 2.0;
	middle = (b.greater_corner.y + b.lesser_corner.y) / 2.0;
	b.greater_corner.y = middle + max / 2.0;
	b.lesser_corner.y = middle - max / 2.0;
	middle = (b.greater_corner.z + b.lesser_corner.z) / 2.0;
	b.greater_corner.z = middle + max / 2.0;
	b.lesser_corner.z = middle - max / 2.0;
	
	//nudge the boundaries up by enough bits so that all particles' keys get generated within the allowed range
	bumpBox(b, HUGE_VAL);
}

} //close namespace SFC

#endif //SFC_H
