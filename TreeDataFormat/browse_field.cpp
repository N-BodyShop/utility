/** @file browse_field.cpp
 Given a field file, browse the values by index.
 */

#include <iostream>
#include <fstream>
#include <vector>

#include "tree_xdr.h"

using namespace std;

void printFieldValue(ostream& os, const FieldHeader& fh, void* data, int i) {
	switch(fh.code) {
		case int8:
			for(unsigned int j = 0; j < fh.dimensions; ++j)
				os << static_cast<char *>(data)[i * fh.dimensions + j] << " ";
			break;
		case uint8:
			for(unsigned int j = 0; j < fh.dimensions; ++j)
				os << static_cast<unsigned char *>(data)[i * fh.dimensions + j] << " ";
			break;
		case int16:
			for(unsigned int j = 0; j < fh.dimensions; ++j)
				os << static_cast<short *>(data)[i * fh.dimensions + j] << " ";
			break;
		case uint16:
			for(unsigned int j = 0; j < fh.dimensions; ++j)
				os << static_cast<unsigned short *>(data)[i * fh.dimensions + j] << " ";
			break;
		case int32:
			for(unsigned int j = 0; j < fh.dimensions; ++j)
				os << static_cast<int *>(data)[i * fh.dimensions + j] << " ";
			break;
		case uint32:
			for(unsigned int j = 0; j < fh.dimensions; ++j)
				os << static_cast<unsigned int *>(data)[i * fh.dimensions + j] << " ";
			break;
		case int64:
			for(unsigned int j = 0; j < fh.dimensions; ++j)
				os << static_cast<int64_t *>(data)[i * fh.dimensions + j] << " ";
			break;
		case uint64:
			for(unsigned int j = 0; j < fh.dimensions; ++j)
				os << static_cast<u_int64_t *>(data)[i * fh.dimensions + j] << " ";
			break;
		case float32:
			for(unsigned int j = 0; j < fh.dimensions; ++j)
				os << static_cast<float *>(data)[i * fh.dimensions + j] << " ";
			break;
		case float64:
			for(unsigned int j = 0; j < fh.dimensions; ++j)
				os << static_cast<double *>(data)[i * fh.dimensions + j] << " ";
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
	
	int i;
	cout << "Which index do you want to see? ";
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
	
	cerr << "Done." << endl;
}
