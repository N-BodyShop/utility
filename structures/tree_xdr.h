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
	static const unsigned int sizeBytes = 76;

	static const int MagicNumber = 3184622;

	int magic;
	
	double time;
	u_int64_t numNodes;
	u_int64_t numParticles;
	OrientedBox<double> boundingBox;
	
	TreeHeader() : magic(MagicNumber) { }
	
	/// Output operator, used for formatted display
	friend std::ostream& operator<< (std::ostream& os, const TreeHeader& h) {
		return os << "Time: " << h.time
				<< "\nTotal number of particles: " << h.numParticles
				<< "\nTotal number of nodes: " << h.numNodes
				<< "\nBounding box: " << h.boundingBox;
	}
};

template <>
inline bool_t xdr_template(XDR* xdrs, TreeHeader* h) {
	return (xdr_template(xdrs, &(h->magic))
			&& xdr_template(xdrs, &(h->time))
			&& xdr_template(xdrs, &(h->numNodes))
			&& xdr_template(xdrs, &(h->numParticles))
			&& xdr_template(xdrs, &(h->boundingBox)));
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

/** Write a field to an XDR stream.  You have to write the header yourself
 before you call this if you want to use the result again.
 */
template <typename T>
bool writeField(XDR* xdrs, const u_int64_t N, T* data) {
	for(u_int64_t i = 0; i < N; ++i) {
		if(!xdr_template(xdrs, data + i)) {
			return false;
		}
	}
	return true;
}

/** Given the type code in the header, write the correct data type
 for the field.
 */
bool writeField(FieldHeader fh, XDR* xdrs, void* data) {
	if(!xdr_template(xdrs, &fh))
		return false;
	switch(fh.code) {
		case int8:
			return writeField(xdrs, fh.dimensions * fh.numParticles, reinterpret_cast<char *>(data));
		case uint8:
			return writeField(xdrs, fh.dimensions * fh.numParticles, reinterpret_cast<unsigned char *>(data));
		case int16:
			return writeField(xdrs, fh.dimensions * fh.numParticles, reinterpret_cast<short *>(data));
		case uint16:
			return writeField(xdrs, fh.dimensions * fh.numParticles, reinterpret_cast<unsigned short *>(data));
		case int32:
			return writeField(xdrs, fh.dimensions * fh.numParticles, reinterpret_cast<int *>(data));
		case uint32:
			return writeField(xdrs, fh.dimensions * fh.numParticles, reinterpret_cast<unsigned int *>(data));
		case int64:
			return writeField(xdrs, fh.dimensions * fh.numParticles, reinterpret_cast<int64_t *>(data));
		case uint64:
			return writeField(xdrs, fh.dimensions * fh.numParticles, reinterpret_cast<u_int64_t *>(data));
		case float32:
			return writeField(xdrs, fh.dimensions * fh.numParticles, reinterpret_cast<float *>(data));
		case float64:
			return writeField(xdrs, fh.dimensions * fh.numParticles, reinterpret_cast<double *>(data));
		default:
			return false;
	}
}

/** Allocate for and read in a field from an XDR stream.
 */
template <typename T>
T* readField(XDR* xdrs, const u_int64_t N) {
	T* data = new T[N];
	for(u_int64_t i = 0; i < N; ++i) {
		if(!xdr_template(xdrs, data + i)) {
			delete[] data;
			return 0;
		}
	}
	return data;
}

/** Given the type code in the header, reads in the correct type of data.
 */
void* readField(const FieldHeader& fh, XDR* xdrs) {
	switch(fh.code) {
		case int8:
			return readField<char>(xdrs, fh.dimensions * fh.numParticles);
		case uint8:
			return readField<unsigned char>(xdrs, fh.dimensions * fh.numParticles);
		case int16:
			return readField<short>(xdrs, fh.dimensions * fh.numParticles);
		case uint16:
			return readField<unsigned short>(xdrs, fh.dimensions * fh.numParticles);
		case int32:
			return readField<int>(xdrs, fh.dimensions * fh.numParticles);
		case uint32:
			return readField<unsigned int>(xdrs, fh.dimensions * fh.numParticles);
		case int64:
			return readField<int64_t>(xdrs, fh.dimensions * fh.numParticles);
		case uint64:
			return readField<u_int64_t>(xdrs, fh.dimensions * fh.numParticles);
		case float32:
			return readField<float>(xdrs, fh.dimensions * fh.numParticles);
		case float64:
			return readField<double>(xdrs, fh.dimensions * fh.numParticles);
		default:
			return 0;
	}
}

// XDR has a minimum size of 4 bytes.
unsigned int mySizeof(DataTypeCode code) {
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
	return xdr_setpos(xdrs, fh.sizeBytes + index * fh.dimensions * mySizeof(fh.code));
}

/** The tree file contains these structures.  With the total number of
 nodes and particles and the bounding box (found in the tree header)
 you can do a tree-based search for particles using these nodes.
 */
struct BasicTreeNode {
	u_int64_t numNodesLeft;
	u_int64_t numParticlesLeft;	
};

template <>
inline bool_t xdr_template(XDR* xdrs, BasicTreeNode* node) {
	return xdr_template(xdrs, &(node->numNodesLeft)) 
			&& xdr_template(xdrs, &(node->numParticlesLeft));
}


#endif //TREE_XDR_H
