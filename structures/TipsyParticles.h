/** @file TipsyParticles.h
 This file defines the particle classes used in Tipsy.
 @author Graeme Lufkin (gwl@u.washington.edu)
 @date Created November 14, 2002
 @version 1.0
 */
 
#ifndef TIPSYPARTICLES_H
#define TIPSYPARTICLES_H

#include <iostream>

#include "Vector3D.h"

namespace Tipsy {

/// The number of dimensions, always 3.
const int MAXDIM = 3;

/** @deprecated A shorthand for an infinite loop. */
#define forever for(;;)

typedef float Real;

class simple_particle {
public:
	static const unsigned int sizeBytes = 28;
	
	/** The mass of the particle. */
    Real mass;
	/** The position of this particle. */
    Vector3D<Real> pos;
	/** The velocity vector of this particle. */
    Vector3D<Real> vel;
	
	simple_particle() : mass(0) { }
	
	bool operator==(const simple_particle& p) const {
		return (mass == p.mass) && (pos == p.pos) && (vel == p.vel);
	}
	
	/// Output operator, used for formatted display
	friend std::ostream& operator<< (std::ostream& os, const simple_particle& p) {
		return os << "Mass: " << p.mass
				<< "\nPosition: " << p.pos
				<< "\nVelocity: " << p.vel;
	}
	
};

/** A particle representing a blob of gas, this is used in SPH calculations. */
class gas_particle : public simple_particle {
public:

	static const unsigned int sizeBytes = 48;
	
	/** The local density of gas at this particle's location. */
    Real rho;
	/** The temperature of the gas in this particle. */
    Real temp;
	/** The gravitational softening length of this particle. */
    Real hsmooth;
	/** The metal content of this gas particle. */
    Real metals;
	/** The gravitational potential at this particle. */
    Real phi;

	/// Default constructor sets all values to zero
	gas_particle() : rho(0), temp(0), hsmooth(0), metals(0), phi(0) { }
	
	bool operator==(const gas_particle& p) const {
		return (mass == p.mass) && (pos == p.pos) && (vel == p.vel) && (rho == p.rho)
				&& (temp == p.temp) && (hsmooth == p.hsmooth) && (metals == p.metals)
				&& (phi == p.phi);
	}
	
	/// Output operator, used for formatted display
	friend std::ostream& operator<< (std::ostream& os, const gas_particle& p) {
		return os << "Mass: " << p.mass
			<< "\nPosition: " << p.pos << "  Magnitude: " << p.pos.length()
			<< "\nVelocity: " << p.vel << "  Magnitude: " << p.vel.length()
			<< "\nRho: " << p.rho
			<< "\nTemperature: " << p.temp
			<< "\nGravitational Softening: " << p.hsmooth
			<< "\nMetals: " << p.metals
			<< "\nPhi: " << p.phi;
	}
};

/** A particle of dark matter, this interacts via gravity only. */
class dark_particle : public simple_particle {
public:
	
	static const unsigned int sizeBytes = 36;
	
	/** The gravitational softening length of this particle. */
    Real eps;
	/** The gravitational potential at this particle. */
    Real phi;

	/// Default constructor sets all values to zero
	dark_particle() : eps(0), phi(0) { }
	
	bool operator==(const dark_particle& p) const {
		return (mass == p.mass) && (pos == p.pos) && (vel == p.vel)
				&& (eps == p.eps) && (phi == p.phi);
	}
	
	/// Output operator, used for formatted display
	friend std::ostream& operator<< (std::ostream& os, const dark_particle& p) {
		return os << "Mass: " << p.mass
			<< "\nPosition: " << p.pos << "  Magnitude: " << p.pos.length()
			<< "\nVelocity: " << p.vel << "  Magnitude: " << p.vel.length()
			<< "\nGravitational Softening: " << p.eps
			<< "\nPhi: " << p.phi;
	}
};

/** A star particle interacts via gravity, and has some extra properties. */
class star_particle : public simple_particle {
public:
	
	static const unsigned int sizeBytes = 44;
	
	/** The metallicity of this star particle. */
    Real metals;
	/** The time of formation of this star particle. */
    Real tform;
	/** The gravitational softening length of this particle. */
    Real eps;
	/** The gravitational potential at this particle. */
    Real phi;

	/// Default constructor sets all values to zero
	star_particle() : metals(0), tform(0), eps(0), phi(0) { }
	
	bool operator==(const star_particle& p) const {
		return (mass == p.mass) && (pos == p.pos) && (vel == p.vel)
				&& (metals == p.metals) && (tform == p.tform) && (eps == p.eps)
				&& (phi == p.phi);
	}

	/// Output operator, used for formatted display
	friend std::ostream& operator<< (std::ostream& os, const star_particle& p) {
		return os << "Mass: " << p.mass
			<< "\nPosition: " << p.pos << "  Magnitude: " << p.pos.length()
			<< "\nVelocity: " << p.vel << "  Magnitude: " << p.vel.length()
			<< "\nMetals: " << p.metals
			<< "\nFormation Time: " << p.tform
			<< "\nGravitational Softening: " << p.eps
			<< "\nPhi: " << p.phi;
	}
};

/** A particle representing a group of other particles.
 This class is an extension, it is not used in Tipsy.
*/
class group_particle {
public:
	/** An ID for this group. */
	int groupNumber;
	/** The total mass of this group. */
	Real total_mass;
	/** The reference vector to a member of the group.  This is present so that periodic boundary conditions can work. */
	Vector3D<Real> reference;
	/** The location, relative to the reference, of the center of mass of this group. */
	Vector3D<Real> cm;
	/** The center of mass velocity of this group. */
	Vector3D<Real> cm_velocity;
	/** The number of particles in this group. */
	int numMembers;
	/** The distance from the center of mass to the farthest member of the group. */
	Real radius;
	
	/// Default constructor sets all values to zero
	group_particle() : groupNumber(0), total_mass(0), numMembers(0), radius(0) { }
	
	bool operator==(const group_particle& p) const {
		return (groupNumber == p.groupNumber) && (total_mass == p.total_mass) && ((cm + reference) == (p.cm + p.reference)) && (cm_velocity == p.cm_velocity)
				&& (numMembers == p.numMembers) && (radius == p.radius);
	}

	/// Output operator, used for formatted display
	friend std::ostream& operator<< (std::ostream& os, const group_particle& p) {
		return os << "Group number: " << p.groupNumber
			<< "\nTotal mass: " << p.total_mass
			<< "\nCenter of mass: " << (p.cm + p.reference)
			<< "\nCenter of mass velocity: " << p.cm_velocity << "  Magnitude: " << p.cm_velocity.length()
			<< "\nNumber of members: " << p.numMembers
			<< "\nRadius to farthest member: " << p.radius;
	}
};		

} //close namespace Tipsy

#endif //TIPSYPARTICLES_H
