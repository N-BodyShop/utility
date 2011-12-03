// Convert a salsa file into Tipsy format.
// Still very incomplete: just reads dark, and doesn't get all the
// fields.

#include <iostream>
#include <cstdio>

#include "tree_xdr.h"
#include "TipsyFile.h"

using namespace std;
using namespace TypeHandling;
using namespace Tipsy;

static double fh_time; // gross, but quick way to get time

// Returns total number of particles of a given type.
// Assumes the position attribute is available
int64_t getCount(char *typedir // directory containing the data for
		 // the type of particle of interest
		 )
{
    const int FILELEN = 256;
    char filename[FILELEN];

    strncpy(filename, typedir, FILELEN);
    strcat(filename, "/position");

    FILE* infile = fopen(filename, "rb");
    if(!infile) {
	cerr << "Couldn't open field file \"" << filename << "\"" << endl;
	return 0;  // Assume there is none of this particle type
	}

    XDR xdrs;
    FieldHeader fh;
    xdrstdio_create(&xdrs, infile, XDR_DECODE);

    if(!xdr_template(&xdrs, &fh)) {
	throw XDRException("Couldn't read header from file!");
	}
    if(fh.magic != FieldHeader::MagicNumber) {
	throw XDRException("This file does not appear to be a field file (magic number doesn't match).");
	}
    if(fh.dimensions != 3) {
	throw XDRException("Wrong dimension of positions.");
	}
    fh_time = fh.time;
    xdr_destroy(&xdrs);
    fclose(infile);
    return fh.numParticles;
    }

void *readFieldData(char *filename, FieldHeader &fh, unsigned int dim)
{
    FILE* infile = fopen(filename, "rb");
    if(!infile) {
	throw XDRException("Couldn't open field file");
	}

    XDR xdrs;
    xdrstdio_create(&xdrs, infile, XDR_DECODE);

    if(!xdr_template(&xdrs, &fh)) {
	throw XDRException("Couldn't read header from file!");
	}
    if(fh.magic != FieldHeader::MagicNumber) {
	throw XDRException("This file does not appear to be a field file (magic number doesn't match).");
	}
    if(fh.dimensions != dim) {
	throw XDRException("Wrong dimension of positions.");
	}

    void* data = readField(fh, &xdrs);
	
    if(data == 0) {
	throw XDRException("Had problems reading in the field");
	}
    xdr_destroy(&xdrs);
    fclose(infile);
    return data;
    }

template <class PartVecT>
void getPos(PartVecT &p, // reference to particle array
	    char *typedir // directory with file
	    )
{
    const int FILELEN = 256;
    char filename[FILELEN];
    FieldHeader fh;

    strncpy(filename, typedir, FILELEN);
    strcat(filename, "/position");
    
    void *data = readFieldData(filename, fh, 3);

    for(unsigned int i = 0; i < fh.numParticles; ++i) {
	switch(fh.code) {
	case float32:
	    for(unsigned int j = 0; j < fh.dimensions; ++j)
		p[i].pos[j] = static_cast<float *>(data)[fh.dimensions * i + j];
	    break;
	case float64:
	    for(unsigned int j = 0; j < fh.dimensions; ++j)
		p[i].pos[j] = static_cast<double *>(data)[fh.dimensions * i + j];
	    break;
	default:
	    throw XDRException("I don't recognize the type of this field!");
	    }
	    }
	
    deleteField(fh, data);
    }

template <class PartVecT>
void getVel(PartVecT &p, // reference to particle array
	    char *typedir // directory with file
	    )
{
    const int FILELEN = 256;
    char filename[FILELEN];
    FieldHeader fh;

    strncpy(filename, typedir, FILELEN);
    strcat(filename, "/velocity");
    
    void *data = readFieldData(filename, fh, 3);

    for(unsigned int i = 0; i < fh.numParticles; ++i) {
	switch(fh.code) {
	case float32:
	    for(unsigned int j = 0; j < fh.dimensions; ++j)
		p[i].vel[j] = static_cast<float *>(data)[fh.dimensions * i + j];
	    break;
	case float64:
	    for(unsigned int j = 0; j < fh.dimensions; ++j)
		p[i].vel[j] = static_cast<double *>(data)[fh.dimensions * i + j];
	    break;
	default:
	    throw XDRException("I don't recognize the type of this field!");
	    }
	    }
	
    deleteField(fh, data);
    }

template <class PartVecT>
void getMass(PartVecT &p, // reference to particle array
	    char *typedir // directory with file
	    )
{
    const int FILELEN = 256;
    char filename[FILELEN];
    FieldHeader fh;

    strncpy(filename, typedir, FILELEN);
    strcat(filename, "/mass");
    
    void *data = readFieldData(filename, fh, 1);

    for(unsigned int i = 0; i < fh.numParticles; ++i) {
	switch(fh.code) {
	case float32:
	    p[i].mass = static_cast<float *>(data)[fh.dimensions * i];
	    break;
	case float64:
	    p[i].mass = static_cast<double *>(data)[fh.dimensions * i];
	    break;
	default:
	    throw XDRException("I don't recognize the type of this field!");
	    }
	    }
	
    deleteField(fh, data);
    }

template <class PartVecT>
void getSoft(PartVecT &p, // reference to particle array
	    char *typedir // directory with file
	    )
{
    const int FILELEN = 256;
    char filename[FILELEN];
    FieldHeader fh;

    strncpy(filename, typedir, FILELEN);
    strcat(filename, "/softening");
    
    void *data = readFieldData(filename, fh, 1);

    for(unsigned int i = 0; i < fh.numParticles; ++i) {
	switch(fh.code) {
	case float32:
	    p[i].eps = static_cast<float *>(data)[fh.dimensions * i];
	    break;
	case float64:
	    p[i].eps = static_cast<double *>(data)[fh.dimensions * i];
	    break;
	default:
	    throw XDRException("I don't recognize the type of this field!");
	    }
	    }
	
    deleteField(fh, data);
    }

template <class PartVecT>
void getPhi(PartVecT &p, // reference to particle array
	    char *typedir // directory with file
	    )
{
    const int FILELEN = 256;
    char filename[FILELEN];
    FieldHeader fh;

    strncpy(filename, typedir, FILELEN);
    strcat(filename, "/potential");
    
    void *data = readFieldData(filename, fh, 1);

    for(unsigned int i = 0; i < fh.numParticles; ++i) {
	switch(fh.code) {
	case float32:
	    p[i].phi = static_cast<float *>(data)[fh.dimensions * i];
	    break;
	case float64:
	    p[i].phi = static_cast<double *>(data)[fh.dimensions * i];
	    break;
	default:
	    throw XDRException("I don't recognize the type of this field!");
	    }
	    }
	
    deleteField(fh, data);
    }

template <class PartVecT>
void getHSmooth(PartVecT &p, // reference to particle array
	    char *typedir // directory with file
	    )
{
    const int FILELEN = 256;
    char filename[FILELEN];
    FieldHeader fh;

    strncpy(filename, typedir, FILELEN);
    strcat(filename, "/softening");
    
    void *data = readFieldData(filename, fh, 1);

    for(unsigned int i = 0; i < fh.numParticles; ++i) {
	switch(fh.code) {
	case float32:
	    p[i].hsmooth = static_cast<float *>(data)[fh.dimensions * i];
	    break;
	case float64:
	    p[i].hsmooth = static_cast<double *>(data)[fh.dimensions * i];
	    break;
	default:
	    throw XDRException("I don't recognize the type of this field!");
	    }
	    }
	
    deleteField(fh, data);
    }

template <class PartVecT>
void getRho(PartVecT &p, // reference to particle array
	    char *typedir // directory with file
	    )
{
    const int FILELEN = 256;
    char filename[FILELEN];
    FieldHeader fh;

    strncpy(filename, typedir, FILELEN);
    strcat(filename, "/density");
    
    void *data = readFieldData(filename, fh, 1);

    for(unsigned int i = 0; i < fh.numParticles; ++i) {
	switch(fh.code) {
	case float32:
	    p[i].rho = static_cast<float *>(data)[fh.dimensions * i];
	    break;
	case float64:
	    p[i].rho = static_cast<double *>(data)[fh.dimensions * i];
	    break;
	default:
	    throw XDRException("I don't recognize the type of this field!");
	    }
	    }
	
    deleteField(fh, data);
    }

template <class PartVecT>
void getTemp(PartVecT &p, // reference to particle array
	    char *typedir // directory with file
	    )
{
    const int FILELEN = 256;
    char filename[FILELEN];
    FieldHeader fh;

    strncpy(filename, typedir, FILELEN);
    strcat(filename, "/temperature");
    
    void *data = readFieldData(filename, fh, 1);

    for(unsigned int i = 0; i < fh.numParticles; ++i) {
	switch(fh.code) {
	case float32:
	    p[i].temp = static_cast<float *>(data)[fh.dimensions * i];
	    break;
	case float64:
	    p[i].temp = static_cast<double *>(data)[fh.dimensions * i];
	    break;
	default:
	    throw XDRException("I don't recognize the type of this field!");
	    }
	    }
	
    deleteField(fh, data);
    }

template <class PartVecT>
void getMetals(PartVecT &p, // reference to particle array
	    char *typedir // directory with file
	    )
{
    const int FILELEN = 256;
    char filename[FILELEN];
    FieldHeader fh;

    strncpy(filename, typedir, FILELEN);
    strcat(filename, "/metals");
    
    void *data = readFieldData(filename, fh, 1);

    for(unsigned int i = 0; i < fh.numParticles; ++i) {
	switch(fh.code) {
	case float32:
	    p[i].metals = static_cast<float *>(data)[fh.dimensions * i];
	    break;
	case float64:
	    p[i].metals = static_cast<double *>(data)[fh.dimensions * i];
	    break;
	default:
	    throw XDRException("I don't recognize the type of this field!");
	    }
	    }
	
    deleteField(fh, data);
    }

template <class PartVecT>
void getTForm(PartVecT &p, // reference to particle array
	    char *typedir // directory with file
	    )
{
    const int FILELEN = 256;
    char filename[FILELEN];
    FieldHeader fh;

    strncpy(filename, typedir, FILELEN);
    strcat(filename, "/formationtime");
    
    void *data = readFieldData(filename, fh, 1);

    for(unsigned int i = 0; i < fh.numParticles; ++i) {
	switch(fh.code) {
	case float32:
	    p[i].tform = static_cast<float *>(data)[fh.dimensions * i];
	    break;
	case float64:
	    p[i].tform = static_cast<double *>(data)[fh.dimensions * i];
	    break;
	default:
	    throw XDRException("I don't recognize the type of this field!");
	    }
	    }
	
    deleteField(fh, data);
    }

int main(int argc, char** argv) {
	if(argc < 3) {
		cerr << "Usage: salsa2tipsy basedir outfile" << endl;
		return 1;
	    }
	FieldHeader fh;
	const int FILELEN = 256;
	char filename[FILELEN];
	
	strncpy(filename, argv[1], FILELEN);
	strcat(filename, "/dark");
	int64_t nDark = getCount(filename);

	strncpy(filename, argv[1], FILELEN);
	strcat(filename, "/gas");
	int64_t nSph = getCount(filename);

	strncpy(filename, argv[1], FILELEN);
	strcat(filename, "/star");
	int64_t nStar= getCount(filename);

	TipsyFile tf("tst", nSph, nDark, nStar);
	
	tf.h.time = fh_time;

	strncpy(filename, argv[1], FILELEN);
	strcat(filename, "/dark");
	getPos(tf.darks, filename);
	getMass(tf.darks, filename);
	getVel(tf.darks, filename);
	getSoft(tf.darks, filename);
	getPhi(tf.darks, filename);

	strncpy(filename, argv[1], FILELEN);
	strcat(filename, "/gas");
	getPos(tf.gas, filename);
	getMass(tf.gas, filename);
	getVel(tf.gas, filename);
	getPhi(tf.gas, filename);
	getHSmooth(tf.gas, filename);
	getRho(tf.gas, filename);
	getTemp(tf.gas, filename);
	getMetals(tf.gas, filename);

	strncpy(filename, argv[1], FILELEN);
	strcat(filename, "/star");
	getPos(tf.stars, filename);
	getMass(tf.stars, filename);
	getVel(tf.stars, filename);
	getSoft(tf.stars, filename);
	getPhi(tf.stars, filename);
	getMetals(tf.stars, filename);
	getTForm(tf.stars, filename);
	
	Tipsy::TipsyWriter w(argv[2], tf.h);
	w.writeHeader();
	w.seekParticleNum(0);
	for(int64_t i = 0; i < nSph; i++)
	    w.putNextGasParticle(tf.gas[i]);
	for(int64_t i = 0; i < nDark; i++)
	    w.putNextDarkParticle(tf.darks[i]);
	for(int64_t i = 0; i < nStar; i++)
	    w.putNextStarParticle(tf.stars[i]);
	    
	cerr << "Done." << endl;	
    }
