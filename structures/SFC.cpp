/** @file SFC.cpp
 Definition of key-making function.
 Must be placed in separate compilation unit because the compiler
 requires an extra flag to prevent the floats from being put into
 registers, which would mess up the calculation, which depends on
 the precise IEEE representation of 32-bit floating point numbers.
 */

#include "SFC.h"

namespace SFC {

/// Out of three floats, make a Morton order (z-order) space-filling curve key
// Cannot make this function inline, or g++ with -O2 or higher generates incorrect code
/** Given the floating point numbers for the location, construct the key. 
 The key uses 21 of the 23 bits for the floats of the x, y, and z coordinates
 of the position vector.  This process will only make sense if the position
 coordinates are in the range [1,2).  The mantissa bits are taken, and interleaved
 in xyz order to form the key.  This makes the key a position on the z-ordering
 space-filling curve. 
 */
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

/** Given a key, create a vector of floats representing a position.
 This is almost the inverse of the makeKey() function.  Since the key
 does not use all the bits, each float generated here will have its last
 two mantissa bits set zero, regardless of the values when the key was
 originally generated.
 */
Vector3D<float> makeVector(Key k) {
	Vector3D<float> v(1.0, 1.0, 1.0);
	int* ix = reinterpret_cast<int *>(&v.x);
	int* iy = reinterpret_cast<int *>(&v.y);
	int* iz = reinterpret_cast<int *>(&v.z);
	for(int mask = (1 << 2); mask < (1 << 23); mask <<= 1) {
		if(k & 4)
			*ix |= mask;
		if(k & 2)
			*iy |= mask;
		if(k & 1)
			*iz |= mask;
		k >>= 3;	
	}
	return v;
}

} //close namespace SFC
