/** @file browse_field.cpp
 Given a field file, browse the values by index.
 */

#include <iostream>
#include <fstream>
#include <vector>

#include "tree_xdr.h"

using namespace std;

template <typename T>
void printFieldValue(ostream& os, T* data, int i) {
	os << data[i];
}

template <typename T>
void printFieldValue(ostream& os, const unsigned int dimensions, T* data, int i) {
	if(dimensions == 3)
		printFieldValue(os, reinterpret_cast<Vector3D<T> *>(data), i);
	else
		printFieldValue(os, data, i);
}

void printFieldValue(ostream& os, const FieldHeader& fh, void* data, int i) {
	switch(fh.code) {
		case int8:
			printFieldValue(os, fh.dimensions, static_cast<char *>(data), i);
			break;
		case uint8:
			printFieldValue(os, fh.dimensions, static_cast<unsigned char *>(data), i);
			break;
		case int16:
			printFieldValue(os, fh.dimensions, static_cast<short *>(data), i);
			break;
		case uint16:
			printFieldValue(os, fh.dimensions, static_cast<unsigned short *>(data), i);
			break;
		case int32:
			printFieldValue(os, fh.dimensions, static_cast<int *>(data), i);
			break;
		case uint32:
			printFieldValue(os, fh.dimensions, static_cast<unsigned int *>(data), i);
			break;
		case int64:
			printFieldValue(os, fh.dimensions, static_cast<int64_t *>(data), i);
			break;
		case uint64:
			printFieldValue(os, fh.dimensions, static_cast<u_int64_t *>(data), i);
			break;
		case float32:
			printFieldValue(os, fh.dimensions, static_cast<float *>(data), i);
			break;
		case float64:
			printFieldValue(os, fh.dimensions, static_cast<double *>(data), i);
			break;
		default:
			os << "I don't recognize the type of this field!";
	}
}

int main(int argc, char** argv) {
	if(argc < 2) {
		cerr << "Usage: browse_field field_file" << endl;
		return 1;
	}
	
	FILE* infile = fopen(argv[1], "rb");
	if(!infile) {
		cerr << "Couldn't open file \"" << argv[1] << "\"" << endl;
		return 2;
	}
	
	XDR xdrs;
	xdrstdio_create(&xdrs, infile, XDR_DECODE);
	FieldHeader fh;
	if(!xdr_template(&xdrs, &fh)) {
		cerr << "Couldn't read header from file!" << endl;
		return 3;
	}
	if(fh.magic != FieldHeader::MagicNumber) {
		cerr << "This file does not appear to be a field file (magic number doesn't match)." << endl;
		return 4;
	}
	
	void* data = readField(fh, &xdrs);
	if(data == 0) {
		cerr << "Had problems reading in the field" << endl;
		return 5;
	}
	
	xdr_destroy(&xdrs);
	fclose(infile);
	
	cout << "File \"" << argv[1] << "\" contains a field.  Header:\n" << fh << endl;
	cout << "Minimum: ";
	printFieldValue(cout, fh, data, fh.numParticles);
	cout << "\nMaximum: ";
	printFieldValue(cout, fh, data, fh.numParticles + 1);
		
	int i;
	cout << "\nWhich index do you want to see? ";
	while((cin >> i) && (i >= 0)) {
		if(i >= fh.numParticles)
			cout << "Index too high, must be between 0 and " << (fh.numParticles - 1);
		else {
			cout << "Value: ";
			printFieldValue(cout, fh, data, i);
		}
		cout << "\nWhich index do you want to see? ";
		cout.flush();
	}
	
	deleteField(fh, data);
	
	cerr << "Done." << endl;
}
