/** @file tree_xdr.h
 Provides class definitions for tree and field headers used in
 tree-based file format.
 */

#ifndef TREE_XDR_H
#define TREE_XDR_H

#include "Vector3D.h"
#include "OrientedBox.h"
#include "xdr_template.h"

/** The header present in every tree file. 
 In conjunction with the positions file, provides all the information 
 needed to perform a tree-based region query on the particle data 
 held in field files.
 */
class TreeHeader {
public:
	static const unsigned int sizeBytes = 84;

	static const int MagicNumber = 3184622;
	
	static const u_int64_t PeriodicFlag = static_cast<u_int64_t>(1) << 0;

	int magic;
	
	double time;
	u_int64_t numNodes;
	u_int64_t numParticles;
	OrientedBox<double> boundingBox;
	u_int64_t flags;
	
	TreeHeader() : magic(MagicNumber), flags(0) { }
	
	/// Output operator, used for formatted display
	friend std::ostream& operator<< (std::ostream& os, const TreeHeader& h) {
		os << "Time: " << h.time
				<< "\nTotal number of particles: " << h.numParticles
				<< "\nTotal number of nodes: " << h.numNodes
				<< "\nBounding box: " << h.boundingBox
				<< "\nFlags: ";
		if(h.flags == 0)
			os << "None";
		else {
			if(h.flags & TreeHeader::PeriodicFlag)
				os << "Periodic ";
		}
		return os;
	}
};

template <>
inline bool_t xdr_template(XDR* xdrs, TreeHeader* h) {
	return (xdr_template(xdrs, &(h->magic))
			&& xdr_template(xdrs, &(h->time))
			&& xdr_template(xdrs, &(h->numNodes))
			&& xdr_template(xdrs, &(h->numParticles))
			&& xdr_template(xdrs, &(h->boundingBox))
			&& xdr_template(xdrs, &(h->flags)));
}

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

/** The header present in every file of attribute data (a field).
 */
class FieldHeader {
public:
	
	static const unsigned int sizeBytes = 28;

	static const int MagicNumber = 1062053;

	int magic;
	double time;
	u_int64_t numParticles;
	unsigned int dimensions; //1 for scalar, 3 for vector
	DataTypeCode code;
	
	FieldHeader() : magic(MagicNumber) { }
	
	/// Output operator, used for formatted display
	friend std::ostream& operator<< (std::ostream& os, const FieldHeader& h) {
		os << "Time: " << h.time
				<< "\nTotal number of particles: " << h.numParticles
				<< "\nDimensions: " << h.dimensions
				<< "\nData Type: ";
		switch(h.code) {
			case int8:
				os << "Signed integer, 8 bits"; break;
			case uint8:
				os << "Unsigned integer, 8 bits"; break;
			case int16:
				os << "Signed integer, 16 bits"; break;
			case uint16:
				os << "Unsigned integer, 16 bits"; break;
			case int32:
				os << "Signed integer, 32 bits"; break;
			case uint32:
				os << "Unsigned integer, 32 bits"; break;
			case int64:
				os << "Signed integer, 64 bits"; break;
			case uint64:
				os << "Unsigned integer, 64 bits"; break;
			case float32:
				os << "Floating point, 32 bits"; break;
			case float64:
				os << "Floating point, 64 bits"; break;
			default:
				os << "Unknown type!"; break;
		}
		return os;
	}
};

template <>
inline bool_t xdr_template(XDR* xdrs, FieldHeader* h) {
	return (xdr_template(xdrs, &(h->magic))
			&& xdr_template(xdrs, &(h->time))
			&& xdr_template(xdrs, &(h->numParticles))
			&& xdr_template(xdrs, &(h->dimensions))
			&& xdr_template(xdrs, reinterpret_cast<enum_t *>(&(h->code))));
}

/** Write a field to an XDR stream.  You have to write the header and min/max 
 yourself before you call this if you want to use the result again.
 */
template <typename T>
inline bool writeField(XDR* xdrs, const u_int64_t N, T* data) {
	for(u_int64_t i = 0; i < N; ++i) {
		if(!xdr_template(xdrs, data + i)) {
			return false;
		}
	}
	return true;
}

/** For three-dimensional types, change cast to a Vector3D of the type. */
template <typename T>
inline bool writeField(XDR* xdrs, const unsigned int dimensions, const u_int64_t N, T* data) {
	if(dimensions == 3)
		return writeField(xdrs, N, reinterpret_cast<Vector3D<T> * >(data));
	else
		return writeField(xdrs, N, data);
}

/** Given a type-code, cast to the appropriate type of pointer. */
inline bool writeFieldSwitch(XDR* xdrs, DataTypeCode code, const unsigned int dimensions, const u_int64_t N, void* data) {
	switch(code) {
		case int8:
			return writeField(xdrs, dimensions, N, reinterpret_cast<char *>(data));
		case uint8:
			return writeField(xdrs, dimensions, N, reinterpret_cast<unsigned char *>(data));
		case int16:
			return writeField(xdrs, dimensions, N, reinterpret_cast<short *>(data));
		case uint16:
			return writeField(xdrs, dimensions, N, reinterpret_cast<unsigned short *>(data));
		case int32:
			return writeField(xdrs, dimensions, N, reinterpret_cast<int *>(data));
		case uint32:
			return writeField(xdrs, dimensions, N, reinterpret_cast<unsigned int *>(data));
		case int64:
			return writeField(xdrs, dimensions, N, reinterpret_cast<int64_t *>(data));
		case uint64:
			return writeField(xdrs, dimensions, N, reinterpret_cast<u_int64_t *>(data));
		case float32:
			return writeField(xdrs, dimensions, N, reinterpret_cast<float *>(data));
		case float64:
			return writeField(xdrs, dimensions, N, reinterpret_cast<double *>(data));
		default:
			return false;
	}
}

/** Given the type code in the header, write the correct data type
 for the field.  You need to pass a correctly filled-in header, the array
 of values, and pointers to the minimum and maximum values.
 */
inline bool writeField(FieldHeader fh, XDR* xdrs, void* data, void* minValue, void* maxValue) {
	if(fh.dimensions != 1 && fh.dimensions != 3)
		return false;
	if(!xdr_template(xdrs, &fh))
		return false;
	return writeFieldSwitch(xdrs, fh.code, fh.dimensions, 1, minValue)
			&& writeFieldSwitch(xdrs, fh.code, fh.dimensions, 1, maxValue)
			&& writeFieldSwitch(xdrs, fh.code, fh.dimensions, fh.numParticles, data);
}

/** Allocate for and read in a field from an XDR stream.  You need to have
 read the header already.  The min/max pair are put at the end of the array.
 */
template <typename T>
inline T* readField(XDR* xdrs, const u_int64_t N) {
	T* data = new T[N + 2];
	//put min/max at the end
	if(!xdr_template(xdrs, data + N) || !xdr_template(xdrs, data + N + 1)) {
		delete[] data;
		return 0;
	}
	if(data[N] == data[N + 1]) { 
		//if all elements are the same, just copy the value into the array
		for(u_int64_t i = 0; i < N; ++i)
			data[i] = data[N];
	} else {
		for(u_int64_t i = 0; i < N; ++i) {
			if(!xdr_template(xdrs, data + i)) {
				delete[] data;
				return 0;
			}
		}
	}
	return data;
}

template <typename T>
inline void* readField(XDR* xdrs, const unsigned int dimensions, const u_int64_t N) {
	if(dimensions == 3)
		return readField<Vector3D<T> >(xdrs, N);
	else
		return readField<T>(xdrs, N);
}

/** Given the type code in the header, reads in the correct type of data.
 */
inline void* readField(const FieldHeader& fh, XDR* xdrs) {
	if(fh.dimensions != 1 && fh.dimensions != 3)
		return 0;
	switch(fh.code) {
		case int8:
			return readField<char>(xdrs, fh.dimensions, fh.numParticles);
		case uint8:
			return readField<unsigned char>(xdrs, fh.dimensions, fh.numParticles);
		case int16:
			return readField<short>(xdrs, fh.dimensions, fh.numParticles);
		case uint16:
			return readField<unsigned short>(xdrs, fh.dimensions, fh.numParticles);
		case int32:
			return readField<int>(xdrs, fh.dimensions, fh.numParticles);
		case uint32:
			return readField<unsigned int>(xdrs, fh.dimensions, fh.numParticles);
		case int64:
			return readField<int64_t>(xdrs, fh.dimensions, fh.numParticles);
		case uint64:
			return readField<u_int64_t>(xdrs, fh.dimensions, fh.numParticles);
		case float32:
			return readField<float>(xdrs, fh.dimensions, fh.numParticles);
		case float64:
			return readField<double>(xdrs, fh.dimensions, fh.numParticles);
		default:
			return 0;
	}
}

/** Using an array of original indices, allocate for and read in a field,
 putting the elements back in the original order. 
 */
template <typename T>
inline T* readField_reorder(XDR* xdrs, const u_int64_t N, const unsigned int* uids) {
	T* data = new T[N + 2];
	//put min/max at the end
	if(!xdr_template(xdrs, data + N) || !xdr_template(xdrs, data + N + 1)) {
		delete[] data;
		return 0;
	}
	if(data[N] == data[N + 1]) { 
		//if all elements are the same, just copy the value into the array
		for(u_int64_t i = 0; i < N; ++i)
			data[i] = data[N];
	} else {
		for(u_int64_t i = 0; i < N; ++i) {
			if(uids[i] < 0 || uids[i] >= N || !xdr_template(xdrs, data + uids[i])) {
				delete[] data;
				return 0;
			}
		}
	}
	return data;
}

template <typename T>
inline void* readField_reorder(XDR* xdrs, const unsigned int dimensions, const u_int64_t N, const unsigned int* uids) {
	if(dimensions == 3)
		return readField_reorder<Vector3D<T> >(xdrs, N, uids);
	else
		return readField_reorder<T>(xdrs, N, uids);
}

/** Using the type code in the header, read in and reorder a field. 
 */
inline void* readField_reorder(const FieldHeader& fh, XDR* xdrs, const unsigned int* uids) {
	if(fh.dimensions != 1 && fh.dimensions != 3)
		return 0;
	switch(fh.code) {
		case int8:
			return readField_reorder<char>(xdrs, fh.dimensions, fh.numParticles, uids);
		case uint8:
			return readField_reorder<unsigned char>(xdrs, fh.dimensions, fh.numParticles, uids);
		case int16:
			return readField_reorder<short>(xdrs, fh.dimensions, fh.numParticles, uids);
		case uint16:
			return readField_reorder<unsigned short>(xdrs, fh.dimensions, fh.numParticles, uids);
		case int32:
			return readField_reorder<int>(xdrs, fh.dimensions, fh.numParticles, uids);
		case uint32:
			return readField_reorder<unsigned int>(xdrs, fh.dimensions, fh.numParticles, uids);
		case int64:
			return readField_reorder<int64_t>(xdrs, fh.dimensions, fh.numParticles, uids);
		case uint64:
			return readField_reorder<u_int64_t>(xdrs, fh.dimensions, fh.numParticles, uids);
		case float32:
			return readField_reorder<float>(xdrs, fh.dimensions, fh.numParticles, uids);
		case float64:
			return readField_reorder<double>(xdrs, fh.dimensions, fh.numParticles, uids);
		default:
			return 0;
	}
}

template <typename T>
inline bool deleteField(const unsigned int dimensions, void* data) {
	if(dimensions == 3)
		delete[] reinterpret_cast<Vector3D<T> *>(data);
	else
		delete[] reinterpret_cast<T *>(data);
	return true;
}

inline bool deleteField(const FieldHeader& fh, void* data) {
	if(fh.dimensions != 1 && fh.dimensions != 3)
		return false;
	switch(fh.code) {
		case int8:
			return deleteField<char>(fh.dimensions, data);
		case uint8:
			return deleteField<unsigned char>(fh.dimensions, data);
		case int16:
			return deleteField<short>(fh.dimensions, data);
		case uint16:
			return deleteField<unsigned short>(fh.dimensions, data);
		case int32:
			return deleteField<int>(fh.dimensions, data);
		case uint32:
			return deleteField<unsigned int>(fh.dimensions, data);
		case int64:
			return deleteField<int64_t>(fh.dimensions, data);
		case uint64:
			return deleteField<u_int64_t>(fh.dimensions, data);
		case float32:
			return deleteField<float>(fh.dimensions, data);
		case float64:
			return deleteField<double>(fh.dimensions, data);
		default:
			return false;
	}
}

// XDR has a minimum size of 4 bytes.
inline unsigned int mySizeof(DataTypeCode code) {
	switch(code) {
		case int8:
		case uint8:
		case int16:
		case uint16:
		case int32:
		case uint32:
			return 4;
		case int64:
		case uint64:
			return 8;
		case float32:
			return 4;
		case float64:
			return 8;
		default:
			return 0;
	}
}

/** Given a header and a stream, seek to the desired location in the field stream.
 */
inline bool_t seekField(const FieldHeader& fh, XDR* xdrs, const u_int64_t index) {
	return xdr_setpos(xdrs, FieldHeader::sizeBytes + (index + 2) * fh.dimensions * mySizeof(fh.code));
}

/** The tree file contains these structures.  With the total number of
 nodes and particles and the bounding box (found in the tree header)
 you can do a tree-based search for particles using these nodes.
 */
struct BasicTreeNode {
	u_int64_t numNodesLeft;
	u_int64_t numParticlesLeft;
	
	bool operator==(const BasicTreeNode& n) const { 
		return (numNodesLeft == n.numNodesLeft) && (numParticlesLeft == n.numParticlesLeft);
	}
};

template <>
inline bool_t xdr_template(XDR* xdrs, BasicTreeNode* node) {
	return xdr_template(xdrs, &(node->numNodesLeft)) 
			&& xdr_template(xdrs, &(node->numParticlesLeft));
}

#ifdef __CHARMC__
#include "pup.h"

inline void operator|(PUP::er& p, TreeHeader& h) {
	p | h.magic;
	p | h.time;
	p | h.numNodes;
	p | h.numParticles;
	p | h.boundingBox;
	p | h.flags;
}

inline void operator|(PUP::er& p, FieldHeader& h) {
	p | h.magic;
	p | h.time;
	p | h.numParticles;
	p | h.dimensions;
	if(p.isUnpacking()) {
		int enum_int;
		p | enum_int;
		h.code = DataTypeCode(enum_int);
	} else
		p | static_cast<int>(h.code);
}

inline void operator|(PUP::er& p, BasicTreeNode& n) {
	p | n.numNodesLeft;
	p | n.numParticlesLeft;
}

#endif //__CHARMC__


#endif //TREE_XDR_H
