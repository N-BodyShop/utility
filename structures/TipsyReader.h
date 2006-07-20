/** @file TipsyReader.h
 A class that reads a tipsy format file from a stream.
 @author Graeme Lufkin (gwl@u.washington.edu)
 @date Created February 12, 2003
 @version 1.0
 */

#ifndef TIPSYREADER_H
#define TIPSYREADER_H

#include <fstream>
#include <string>
#include <vector>

#include "TipsyParticles.h"
#include "Vector3D.h"
#include "xdr_template.h"

namespace Tipsy {

/** The header in a tipsy format file. */
class header {
public:	
	static const unsigned int sizeBytes = 28; // Only for an x86 machine
	
	/// The time of the output
    double time;
	/// The number of particles of all types in this file
    int nbodies;
	/// The number of dimensions, must be equal to MAXDIM
    int ndim;
	/// The number of SPH (gas) particles in this file
    int nsph;
	/// The number of dark matter particles in this file
    int ndark;
	/// The number of star particles in this file
    int nstar;
    //int pad; //unused on x86
	
	header(int nGas = 0, int nDark = 0, int nStar = 0) : time(0), nbodies(nGas + nDark + nStar), ndim(MAXDIM), nsph(nGas), ndark(nDark), nstar(nStar) { }
	
	/// Output operator, used for formatted display
	friend std::ostream& operator<<(std::ostream& os, const header& h) {
		return os << "Time: " << h.time
			<< "\nnBodies: " << h.nbodies
			<< "\nnDim: " << h.ndim
			<< "\nnSPH: " << h.nsph
			<< "\nnDark: " << h.ndark
			<< "\nnStar: " << h.nstar;
	}
};

/** A particle-at-a-time reader for Tipsy files.
 On little-endian machines, can read native and standard format files.
 On big-endian machine, can read standard format files. 
 */
class TipsyReader {
	bool native;
	bool ok;
	
	int numGasRead, numDarksRead, numStarsRead;
	
	header h;
	
	bool responsible;
 public:
	std::istream* tipsyStream;
	
	bool loadHeader();
	
	//copy constructor and equals operator are private to prevent their use (use takeOverStream instead)
	TipsyReader(const TipsyReader& r);
	TipsyReader& operator=(const TipsyReader& r);
	
public:
	
	TipsyReader() : native(true), ok(false), responsible(false), tipsyStream(0) { }

	/// Load from a file
	TipsyReader(const std::string& filename) : ok(false), responsible(true) {
		tipsyStream = new std::ifstream(filename.c_str(), std::ios::in | std::ios::binary);
		loadHeader();
	}
	
	/** Load from a stream.
	 @note When using the stream interface, the given stream's buffer must
	 outlast your use of the reader object, since the original stream owns
	 the buffer.  For the standard stream cin this is guaranteed.
	 */
	TipsyReader(std::istream& is) : ok(false), responsible(true) {
		tipsyStream = new std::istream(is.rdbuf());
		loadHeader();
	}
	
	/// Use this instead of a copy constructor
	void takeOverStream(TipsyReader& r) {
		native = r.native;
		ok = r.ok;
		numGasRead = r.numGasRead;
		numDarksRead = r.numDarksRead;
		numStarsRead = r.numStarsRead;
		h = r.h;
		if(responsible)
			delete tipsyStream;
		responsible = true;
		r.responsible = false;
		tipsyStream = r.tipsyStream;
	}
	
	/** Reload from a file.
	 */
	bool reload(const std::string& filename) {
		if(responsible)
			delete tipsyStream;
		tipsyStream = new std::ifstream(filename.c_str(), std::ios::in | std::ios::binary);
		responsible = true;
		return loadHeader();
	}
	
	/** Reload from a stream.
	 @note When using the stream interface, the given stream's buffer must
	 outlast your use of the reader object, since the original stream owns
	 the buffer.  For the standard stream cin this is guaranteed.
	 */
	bool reload(std::istream& is) {
		if(responsible)
			delete tipsyStream;
		tipsyStream = new std::istream(is.rdbuf());
		responsible = true;
		return loadHeader();		
	}
	
	~TipsyReader() {
		if(responsible)
			delete tipsyStream;
	}
	
	header getHeader() {
		return h;
	}
	
	bool getNextSimpleParticle(simple_particle& p);
	
	bool getNextGasParticle(gas_particle& p);
	bool getNextDarkParticle(dark_particle& p);
	bool getNextStarParticle(star_particle& p);
	
	bool readAllParticles(gas_particle* gas, dark_particle* darks, star_particle* stars);
	bool readAllParticles(std::vector<gas_particle>& gas, std::vector<dark_particle>& darks, std::vector<star_particle>& stars);
	
	/// Is this file in native byte-order
	bool isNative() const {
		return native;
	}
	
	bool status() const {
		return ok;
	}
	
	bool seekParticleNum(unsigned int num);
	
	bool skipParticles(unsigned int num);
	
	/// Returns the index of the next particle to be read.
	unsigned int tellParticleNum() const {
		return numGasRead + numDarksRead + numStarsRead;
	}
};

class TipsyWriter {
	bool native;
	bool ok;
	FILE *tipsyFp;
	XDR xdrs;
	header h;
	
 public:
	
	bool writeHeader();
	
	/// Write to a file
	TipsyWriter(const std::string& filename, header &parh,
		    bool pnative=false)
	    : native(pnative), ok(false), h(parh) {
	    tipsyFp = fopen(filename.c_str(), "a");  // Create file
	    if(tipsyFp == NULL) {
		ok = false;
		assert(0);
		return;
		}
	    fclose(tipsyFp);
	    tipsyFp = fopen(filename.c_str(), "rb+");
	    if(tipsyFp == NULL) {
		ok = false;
		assert(0);
		return;
		}
	    if(!native)
		xdrstdio_create(&xdrs, tipsyFp, XDR_ENCODE);
	    ok = true;
	}
	
	~TipsyWriter() {
	    if(!native)
		xdr_destroy(&xdrs);
	    fclose(tipsyFp);
	}
	
	bool putNextGasParticle(gas_particle& p);
	bool putNextDarkParticle(dark_particle& p);
	bool putNextStarParticle(star_particle& p);
	
	//	bool writeAllParticles(gas_particle* gas, dark_particle* darks, star_particle* stars);
	//	bool writeAllParticles(std::vector<gas_particle>& gas, std::vector<dark_particle>& darks, std::vector<star_particle>& stars);
	
	/// Is this file in native byte-order
	bool isNative() const {
		return native;
	}
	
	bool status() const {
		return ok;
	}
	
	bool seekParticleNum(unsigned int num);
};

} //close namespace Tipsy

#endif //TIPSYREADER_H
