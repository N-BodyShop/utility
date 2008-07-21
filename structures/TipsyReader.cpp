/** @file TipsyReader.cpp
 @author Graeme Lufkin (gwl@u.washington.edu)
 @date Created February 12, 2003
 @version 1.0
 */

#include <assert.h>
#include "config.h"
#include "xdr_template.h"

#include "TipsyReader.h"

/// XDR conversion for the header structure
inline bool_t xdr_template(XDR* xdrs, Tipsy::header* val) {
	return (xdr_template(xdrs, &(val->time))
			&& xdr_template(xdrs, &(val->nbodies))
			&& xdr_template(xdrs, &(val->ndim)) 
			&& xdr_template(xdrs, &(val->nsph))
			&& xdr_template(xdrs, &(val->ndark))
			&& xdr_template(xdrs, &(val->nstar)));
}

///XDR conversions for the particle types
inline bool_t xdr_template(XDR* xdrs, Tipsy::simple_particle* p) {
	return (xdr_template(xdrs, &(p->mass))
		&& xdr_template(xdrs, &(p->pos))
		&& xdr_template(xdrs, &(p->vel)));
}

inline bool_t xdr_template(XDR* xdrs, Tipsy::gas_particle* p) {
	return (xdr_template(xdrs, &(p->mass))
		&& xdr_template(xdrs, &(p->pos))
		&& xdr_template(xdrs, &(p->vel))
		&& xdr_template(xdrs, &(p->rho))
		&& xdr_template(xdrs, &(p->temp))
		&& xdr_template(xdrs, &(p->hsmooth))
		&& xdr_template(xdrs, &(p->metals))
		&& xdr_template(xdrs, &(p->phi)));
}

inline bool_t xdr_template(XDR* xdrs, Tipsy::dark_particle* p) {
	return (xdr_template(xdrs, &(p->mass))
		&& xdr_template(xdrs, &(p->pos))
		&& xdr_template(xdrs, &(p->vel))
		&& xdr_template(xdrs, &(p->eps))
		&& xdr_template(xdrs, &(p->phi)));
}

inline bool_t xdr_template(XDR* xdrs, Tipsy::star_particle* p) {
	return (xdr_template(xdrs, &(p->mass))
		&& xdr_template(xdrs, &(p->pos))
		&& xdr_template(xdrs, &(p->vel))
		&& xdr_template(xdrs, &(p->metals))
		&& xdr_template(xdrs, &(p->tform))
		&& xdr_template(xdrs, &(p->eps))
		&& xdr_template(xdrs, &(p->phi)));
}

namespace Tipsy {

bool TipsyReader::loadHeader() {
	ok = false;
	
	if(!(*tipsyStream))
		return false;
	
	//read the header in
	tipsyStream->read(reinterpret_cast<char *>(&h), sizeof(h));
	if(!*tipsyStream)
		return false;
	
	if(h.ndim != MAXDIM) { //perhaps it's XDR
		XDR xdrs;
		xdrmem_create(&xdrs, reinterpret_cast<char *>(&h), header::sizeBytes, XDR_DECODE);
		if(!xdr_template(&xdrs, &h) || h.ndim != MAXDIM) { //wasn't xdr format either			
			h.nbodies = h.nsph = h.ndark = h.nstar = 0;
			xdr_destroy(&xdrs);
			return false;
		}
		xdr_destroy(&xdrs);
		
		native = false;
		// xdr format has an integer pad in the header,
		// which we don't need, but must skip.
		// However, if the native format has the pad (the compiler
		// pads structures, then we already took in the pad above.
		if(sizeof(h) != 32) {
		    int pad = 0;
                    assert(sizeof(h) == 28); // Assume that native format
					     // is unpadded.
		    tipsyStream->read(reinterpret_cast<char *>(&pad), 4);
		    }
		
		if(!*tipsyStream)
			return false;
		
	} else
		native = true;
	
	numGasRead = numDarksRead = numStarsRead = 0;
	ok = true;
	return ok;
}

/** Read the next particle, and get the simple_particle representation of it.
 Returns true if successful.  Returns false if the read failed, or no more particles.
 */
bool TipsyReader::getNextSimpleParticle(simple_particle& p) {
	if(!ok)
		return false;
	
	if(numGasRead < h.nsph) {
		++numGasRead;
		gas_particle gp;
		tipsyStream->read(reinterpret_cast<char *>(&gp), gas_particle::sizeBytes);
		if(!(*tipsyStream))
			return false;
		p = gp;
	} else if(numDarksRead < h.ndark) {
		++numDarksRead;
		dark_particle dp;
		tipsyStream->read(reinterpret_cast<char *>(&dp), dark_particle::sizeBytes);
		if(!(*tipsyStream))
			return false;
		p = dp;
	} else if(numStarsRead < h.nstar) {
		++numStarsRead;
		star_particle sp;
		tipsyStream->read(reinterpret_cast<char *>(&sp), star_particle::sizeBytes);
		if(!(*tipsyStream))
			return false;
		p = sp;
	} else //already read all the particles!
		return false;
	
	if(!native) {
		XDR xdrs;
		xdrmem_create(&xdrs, reinterpret_cast<char *>(&p), simple_particle::sizeBytes, XDR_DECODE);
		if(!xdr_template(&xdrs, &p))
			return false;
		xdr_destroy(&xdrs);
	}
	
	return true;
}

/** Get the next gas particle.
 Returns false if the read failed, or already read all the gas particles in this file.
 */
bool TipsyReader::getNextGasParticle(gas_particle& p) {
	if(!ok || !(*tipsyStream))
		return false;
	
	if(numGasRead < h.nsph) {
		++numGasRead;
		tipsyStream->read(reinterpret_cast<char *>(&p), gas_particle::sizeBytes);
		if(numGasRead < h.nsph && !(*tipsyStream))
			return false;
		if(!native) {
			XDR xdrs;
			xdrmem_create(&xdrs, reinterpret_cast<char *>(&p), gas_particle::sizeBytes, XDR_DECODE);
			if(!xdr_template(&xdrs, &p))
				return false;
			xdr_destroy(&xdrs);
		}
	} else
		return false;
	
	return true;
}

/** Get the next dark particle.
 Returns false if the read failed, or already read all the dark particles in this file.
 */
bool TipsyReader::getNextDarkParticle(dark_particle& p) {
	if(!ok || !(*tipsyStream))
		return false;
	
	if(numGasRead != h.nsph) { //some gas still not read, skip them
		if(!seekParticleNum(h.nsph))
			return false;
		numGasRead = h.nsph;
	}
	
	if(numDarksRead < h.ndark) {
		++numDarksRead;
		tipsyStream->read(reinterpret_cast<char *>(&p), dark_particle::sizeBytes);
	// Hack to fix end of stream problem on Macs --trq
		if(numDarksRead < h.ndark && !(*tipsyStream))
			return false;
		if(!native) {
			XDR xdrs;
			xdrmem_create(&xdrs, reinterpret_cast<char *>(&p), dark_particle::sizeBytes, XDR_DECODE);
			if(!xdr_template(&xdrs, &p))
				return false;
			xdr_destroy(&xdrs);
		}
	} else
		return false;
	
	return true;
}

/** Get the next star particle.
 Returns false if the read failed, or already read all the star particles in this file.
 */
bool TipsyReader::getNextStarParticle(star_particle& p) {
	if(!ok || !(*tipsyStream))
		return false;
	
	if(numGasRead != h.nsph || numDarksRead != h.ndark) { //some gas and dark still not read, skip them
		if(!seekParticleNum(h.nsph + h.ndark))
			return false;
		numGasRead = h.nsph;
		numDarksRead = h.ndark;
	}
	
	if(numStarsRead < h.nstar) {
		++numStarsRead;
		tipsyStream->read(reinterpret_cast<char *>(&p), star_particle::sizeBytes);
	// Hack to fix end of stream problem on Macs --trq
		if(numStarsRead < h.nstar && !(*tipsyStream))
			return false;
		if(!native) {
			XDR xdrs;
			xdrmem_create(&xdrs, reinterpret_cast<char *>(&p), star_particle::sizeBytes, XDR_DECODE);
			if(!xdr_template(&xdrs, &p))
				return false;
			xdr_destroy(&xdrs);
		}
	} else
		return false;
	
	return true;
}

/** Read all the particles into arrays.
 Returns false if the read fails.
 */
bool TipsyReader::readAllParticles(gas_particle* gas, dark_particle* darks, star_particle* stars) {
	if(!ok || !(*tipsyStream))
		return false;
	
	//go back to the beginning of the file
	if(!seekParticleNum(0))
		return false;
	
	gas = new gas_particle[h.nsph];
	for(int i = 0; i < h.nsph; ++i) {
		if(!getNextGasParticle(gas[i]))
			return false;
	}
	darks = new dark_particle[h.ndark];
	for(int i = 0; i < h.ndark; ++i) {
		if(!getNextDarkParticle(darks[i]))
			return false;
	}
	stars = new star_particle[h.nstar];
	for(int i = 0; i < h.nstar; ++i) {
		if(!getNextStarParticle(stars[i]))
			return false;
	}
	
	return true;
}

/** Read all the particles into vectors.
 Returns false if the read fails.
 */
bool TipsyReader::readAllParticles(std::vector<gas_particle>& gas, std::vector<dark_particle>& darks, std::vector<star_particle>& stars) {
	if(!ok || !(*tipsyStream))
		return false;
	
	//go back to the beginning of the file
	if(!seekParticleNum(0))
		return false;
	
	gas.resize(h.nsph);
	for(int i = 0; i < h.nsph; ++i) {
		if(!getNextGasParticle(gas[i]))
			return false;
	}
	darks.resize(h.ndark);
	for(int i = 0; i < h.ndark; ++i) {
		if(!getNextDarkParticle(darks[i]))
			return false;
	}
	stars.resize(h.nstar);
	for(int i = 0; i < h.nstar; ++i) {
		if(!getNextStarParticle(stars[i]))
			return false;
	}
	
	return true;
}

/** Seek to the requested particle in the file.
 Returns false if the file seek fails, or if you request too large a particle.
 */
bool TipsyReader::seekParticleNum(unsigned int num) {
	int padSize = 0;
	if(!native)
		padSize = 4;
	else
		padSize = sizeof(h) - header::sizeBytes;

	unsigned int preface = header::sizeBytes + padSize;
	std::streampos seek_position;
	if(num < h.nsph) {
		seek_position = preface + num * (std::streampos) gas_particle::sizeBytes;
		tipsyStream->seekg(seek_position);
		if(!(*tipsyStream))
			return false;
		numGasRead = num;
		numDarksRead = 0;
		numStarsRead = 0;
	} else if(num < (h.nsph + h.ndark)) {
		seek_position = preface + h.nsph * (std::streampos) gas_particle::sizeBytes + (num - h.nsph) * (std::streampos) dark_particle::sizeBytes;
		tipsyStream->seekg(seek_position);
		if(!(*tipsyStream))
			return false;
		numGasRead = h.nsph;
		numDarksRead = num - h.nsph;
		numStarsRead = 0;
	} else if(num < (h.nsph + h.ndark + h.nstar)) {
		seek_position = preface + h.nsph * (std::streampos) gas_particle::sizeBytes + h.ndark * (std::streampos) dark_particle::sizeBytes + (num - h.ndark - h.nsph) * (std::streampos) star_particle::sizeBytes;
		tipsyStream->seekg(seek_position);
		if(!(*tipsyStream))
			return false;
		numGasRead = h.nsph;
		numDarksRead = h.ndark;
		numStarsRead = num - h.ndark - h.nsph;
	} else
		return false;
	
	return true;
}

bool TipsyReader::skipParticles(unsigned int num) {
	for(unsigned int skipped = 0; skipped < num; ++skipped) {
		if(numGasRead < h.nsph) {
			++numGasRead;
			gas_particle gp;
			tipsyStream->read(reinterpret_cast<char *>(&gp),  gas_particle::sizeBytes);
			if(!(*tipsyStream))
				return false;
		} else if(numDarksRead < h.ndark) {
			++numDarksRead;
			dark_particle dp;
			tipsyStream->read(reinterpret_cast<char *>(&dp),  dark_particle::sizeBytes);
			if(!(*tipsyStream))
				return false;
		} else if(numStarsRead < h.nstar) {
			++numStarsRead;
			star_particle sp;
			tipsyStream->read(reinterpret_cast<char *>(&sp),  star_particle::sizeBytes);
			if(!(*tipsyStream))
				return false;
		} else
			return false;
	}
	return true;
}

bool TipsyWriter::writeHeader() {
	
	if(!tipsyFp)
	    return false;
	
	if(native) {
	    fwrite(&h, sizeof(h), 1, tipsyFp);
	    }
	else {
	    xdr_template(&xdrs, &h);
	    }
	return true;
}

/** Write the next gas particle.
 Returns false if the write failed.
 */
bool TipsyWriter::putNextGasParticle(gas_particle& p) {
	if(!tipsyFp)
		return false;
	
	if(native) {
	    fwrite(&p, gas_particle::sizeBytes, 1, tipsyFp);
	    }
	else {
	    if(!xdr_template(&xdrs, &p))
		return false;
	    }
	
	return true;
}

bool TipsyWriter::putNextDarkParticle(dark_particle& p) {
	if(!tipsyFp)
		return false;
	
	if(native) {
	    fwrite(&p, dark_particle::sizeBytes, 1, tipsyFp);
	    }
	else {
	    if(!xdr_template(&xdrs, &p))
		return false;
	    }
	
	return true;
}

bool TipsyWriter::putNextStarParticle(star_particle& p) {
	if(!tipsyFp)
		return false;
	
	if(native) {
	    fwrite(&p, star_particle::sizeBytes, 1, tipsyFp);
	    }
	else {
	    if(!xdr_template(&xdrs, &p))
		return false;
	    }
	
	return true;
}

/** Seek to the requested particle in the file.
 Returns false if the file seek fails, or if you request too large a particle.
 */
bool TipsyWriter::seekParticleNum(unsigned int num) {
	int padSize = 0;
	if(!native)
		padSize = 4;
	else
		padSize = sizeof(h) - header::sizeBytes;

	unsigned int preface = header::sizeBytes + padSize;
	int64_t seek_position;
	xdr_destroy(&xdrs);
	
	if(num < (unsigned int) h.nsph) {
		seek_position = preface + num * gas_particle::sizeBytes;
	} else if(num < (h.nsph + h.ndark)) {
		seek_position = preface + h.nsph * gas_particle::sizeBytes + (num - h.nsph) * dark_particle::sizeBytes;
	} else if(num < (h.nsph + h.ndark + h.nstar)) {
		seek_position = preface + h.nsph * gas_particle::sizeBytes + h.ndark * dark_particle::sizeBytes + (num - h.ndark - h.nsph) * star_particle::sizeBytes;
	} else
		return false;
	fseek(tipsyFp, seek_position, 0);
	xdrstdio_create(&xdrs, tipsyFp, XDR_ENCODE);
	
	return true;
}
} //close namespace Tipsy
