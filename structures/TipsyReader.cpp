/** @file TipsyReader.cpp
 @author Graeme Lufkin (gwl@u.washington.edu)
 @date Created February 12, 2003
 @version 1.0
 */

#include <fstream>

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
	tipsyStream->read(reinterpret_cast<char *>(&h), header::sizeBytes);
	if(!*tipsyStream)
		return false;
	
	if(h.ndim != MAXDIM) { //perhaps it's XDR
		//endian-ness check
		unsigned int bob = 3;
		unsigned char* c = reinterpret_cast<unsigned char *>(&bob);
		if(c[3] == bob) //we're big-endian, can't do xdr from little-endian
			return false;
		
		XDR xdrs;
		xdrmem_create(&xdrs, reinterpret_cast<char *>(&h), header::sizeBytes, XDR_DECODE);
		if(!xdr_template(&xdrs, &h) || h.ndim != MAXDIM) { //wasn't xdr format either			
			h.nbodies = h.nsph = h.ndark = h.nstar = 0;
			return false;
		}
		native = false;
		//xdr format has an integer pad in the header, which we don't need, but must skip
		int pad = 0;
		tipsyStream->read(reinterpret_cast<char *>(&pad), 4);
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
		if(!(*tipsyStream))
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
		if(!(*tipsyStream))
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
		if(!(*tipsyStream))
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
	//endian-ness check
	unsigned int bob = 3;
	unsigned char* c = reinterpret_cast<unsigned char *>(&bob);
	if(c[3] == bob || !native) //we're big-endian, can't do xdr from little-endian
		padSize = 4;

	unsigned int preface = header::sizeBytes + padSize;
	if(num < h.nsph) {
		tipsyStream->seekg(preface + num * gas_particle::sizeBytes);
		if(!(*tipsyStream))
			return false;
		numGasRead = num;
		numDarksRead = 0;
		numStarsRead = 0;
	} else if(num < h.nsph + h.ndark) {
		tipsyStream->seekg(preface + h.nsph * gas_particle::sizeBytes + (num - h.nsph) * dark_particle::sizeBytes);
		if(!(*tipsyStream))
			return false;
		numGasRead = h.nsph;
		numDarksRead = num - h.nsph;
		numStarsRead = 0;
	} else if(num < h.nsph + h.ndark + h.nstar) {
		tipsyStream->seekg(preface + h.nsph * gas_particle::sizeBytes + h.ndark * dark_particle::sizeBytes + (num - h.ndark - h.nsph) * star_particle::sizeBytes);
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
	for(int skipped = 0; skipped < num; ++skipped) {
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

} //close namespace Tipsy
