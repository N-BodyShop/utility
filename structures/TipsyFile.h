/** \file TipsyFile.h
 \author Graeme Lufkin (gwl@u.washington.edu)
 \date Created April 30, 2002
 \version 2.0
 \todo Fully document this file
 */
 
#ifndef TIPSYFILE_H
#define TIPSYFILE_H

#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

#include <rpc/xdr.h>

#include "Vector3D.h"
#include "OrientedBox.h"

using std::string;
using std::cerr;
using std::endl;

/// The number of dimensions, always 3.
const int MAXDIM = 3;

/** \deprecated A shorthand for an infinite loop. */
#define forever for(;;)

typedef float Real;

/** A particle representing a blob of gas, this is used in SPH calculations. */
class gas_particle {
public:

	/** The mass of the gas particle. */
    Real mass;
	/** The position of this particle. */
    Vector3D<Real> pos;
	/** The velocity vector of this particle. */
    Vector3D<Real> vel;
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
	gas_particle() : mass(0), pos(), vel(), rho(0), temp(0), hsmooth(0), metals(0), phi(0) { }
	
	bool operator==(const gas_particle& p) {
		return (mass == p.mass) && (pos == p.pos) && (vel == p.vel) && (rho == p.rho)
				&& (temp == p.temp) && (hsmooth == p.hsmooth) && (metals == p.metals)
				&& (phi == p.phi);
	}
	
	/// Output operator, used for formatted display
	friend std::ostream& operator<< (std::ostream& os, const gas_particle& p) {
		os << "Mass: " << p.mass;
		os << "\nPosition: " << p.pos << "  Magnitude: " << p.pos.length();
		os << "\nVelocity: " << p.vel << "  Magnitude: " << p.vel.length();
		os << "\nRho: " << p.rho;
		os << "\nTemperature: " << p.temp;
		os << "\nGravitational Softening: " << p.hsmooth;
		os << "\nMetals: " << p.metals;
		os << "\nPhi: " << p.phi;
		return os << "\n";
	}
};

typedef gas_particle* gas_particles;

/** A particle of dark matter, this interacts via gravity only. */
class dark_particle {
public:
	
	/** The mass of this dark matter particle. */
    Real mass;
	/** The position of this particle. */
    Vector3D<Real> pos;
	/** The velocity vector of this particle. */
    Vector3D<Real> vel;
	/** The gravitational softening length of this particle. */
    Real eps;
	/** The gravitational potential at this particle. */
    Real phi;

	/// Default constructor sets all values to zero
	dark_particle() : mass(0), pos(), vel(), eps(0), phi(0) { }
	
	bool operator==(const dark_particle& p) {
		return (mass == p.mass) && (pos == p.pos) && (vel == p.vel)
				&& (eps == p.eps) && (phi == p.phi);
	}
	
	/// Output operator, used for formatted display
	friend std::ostream& operator<< (std::ostream& os, const dark_particle& p) {
		os << "Mass: " << p.mass;
		os << "\nPosition: " << p.pos << "  Magnitude: " << p.pos.length();
		os << "\nVelocity: " << p.vel << "  Magnitude: " << p.vel.length();
		os << "\nGravitational Softening: " << p.eps;
		os << "\nPhi: " << p.phi;
		return os << "\n";
	}
};

typedef dark_particle* dark_particles;

class star_particle {
public:
	
	/** The mass of this star particle. */
    Real mass;
	/** The position of this particle. */
    Vector3D<Real> pos;
	/** The velocity vector of this particle. */
    Vector3D<Real> vel;
	/** The metallicity of this star particle. */
    Real metals;
	/** The time of formation of this star particle. */
    Real tform;
	/** The gravitational softening length of this particle. */
    Real eps;
	/** The gravitational potential at this particle. */
    Real phi;

	/// Default constructor sets all values to zero
	star_particle() : mass(0), pos(), vel(), metals(0), tform(0), eps(0), phi(0) { }
	
	bool operator==(const star_particle& p) {
		return (mass == p.mass) && (pos == p.pos) && (vel == p.vel)
				&& (metals == p.metals) && (tform == p.tform) && (eps == p.eps)
				&& (phi == p.phi);
	}

	/// Output operator, used for formatted display
	friend std::ostream& operator<< (std::ostream& os, const star_particle& p) {
		os << "Mass: " << p.mass;
		os << "\nPosition: " << p.pos << "  Magnitude: " << p.pos.length();
		os << "\nVelocity: " << p.vel << "  Magnitude: " << p.vel.length();
		os << "\nMetals: " << p.metals;
		os << "\nFormation Time: " << p.tform;
		os << "\nGravitational Softening: " << p.eps;
		os << "\nPhi: " << p.phi;
		return os << "\n";
	}
};

typedef star_particle* star_particles;

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
	group_particle() : groupNumber(0), total_mass(0), reference(), cm(), cm_velocity(), numMembers(0), radius(0) { }
	
	bool operator==(const group_particle& p) {
		return (groupNumber == p.groupNumber) && (total_mass == p.total_mass) && ((cm + reference) == (p.cm + p.reference)) && (cm_velocity == p.cm_velocity)
				&& (numMembers == p.numMembers) && (radius == p.radius);
	}

	/// Output operator, used for formatted display
	friend std::ostream& operator<< (std::ostream& os, const group_particle& p) {
		os << "Group number: " << p.groupNumber;
		os << "\nTotal mass: " << p.total_mass;
		os << "\nCenter of mass: " << (p.cm + p.reference);
		os << "\nCenter of mass velocity: " << p.cm_velocity << "  Magnitude: " << p.cm_velocity.length();
		os << "\nNumber of members: " << p.numMembers;
		os << "\nRadius to farthest member: " << p.radius;
		return os << "\n";
	}
};
	

/** The header in a tipsy format file. */
struct header {
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
};

/** This class represents a tipsy format file in memory. */
class TipsyFile {
private:

	/// Load the file from disk
	void loadfile();

protected:
	//the xdr stream handler
	XDR xdrs;

	/// XDR function that changes a header structure
	int xdr_convert_header(header &h);
	/// XDR function that changes a gas particle
	int xdr_convert_gas(gas_particle &g);
	/// XDR function that changes a dark matter particle
	int xdr_convert_dark(dark_particle &d);
	/// XDR function that changes a star particle
	int xdr_convert_star(star_particle &s);
	
	bool native;

	bool success;
	
public:
	/// The filename of this tipsy file
	string filename;

	/// The header of this tipsy file
	header h;

	/// The array of gas particles
	gas_particle* gas;
	/// The array of dark matter particles
	dark_particle* darks;
	/// The array of star particles
	star_particle* stars;
		
	TipsyFile() : native(true), success(false), gas(NULL), darks(NULL), stars(NULL) { }

	/// Create a blank tipsy file with specified number of particles
	TipsyFile(const string& fn, int nGas, int nDark = 0, int nStar = 0);
	
	/// Load from a file
	TipsyFile(const string& fn);

	/// Copy constructor
	TipsyFile(const TipsyFile& tf) { 
		h = tf.h;
		native = tf.native;
		success = tf.success;
		filename = tf.filename;
		gas = tf.gas;
		darks = tf.darks;
		stars = tf.stars;
	}
	
	/// Deallocate the memory we used in this file
	~TipsyFile();

	/// Save the current state to disk
	void save() const;
	
	/// Is this file stored in native byte-order?
	bool isNative() const { return native; }

	/// Did the file load successfully?
	bool loadedSuccessfully() const { return success; }

};

class PartialTipsyFile : public TipsyFile {
private:
		
	void loadPartial();

public:
	
	//the header for the full file is stored here.  The header h is the 
	//header structure for the partial file
	header fullHeader;

	int beginParticle;
	int endParticle;
	
	PartialTipsyFile() { }
	PartialTipsyFile(const string& fn, int begin, int end);
	
	void save() const { }
};

class TipsyStats {
private:
	TipsyFile* tf;
public:
	TipsyStats() { }
	TipsyStats(TipsyFile* tf);
	~TipsyStats() { }
	
	double total_mass;
	double gas_mass, dark_mass, star_mass;
	double volume;
	double density;
	OrientedBox<double> bounding_box;
	
	Vector center;
	Vector size;
	
	Vector center_of_mass;
	Vector gas_com;
	Vector dark_com;
	Vector star_com;
	Vector center_of_mass_velocity;
	Vector gas_com_vel;
	Vector dark_com_vel;
	Vector star_com_vel;
	Vector angular_momentum;
	Vector gas_ang_mom;
	Vector dark_ang_mom;
	Vector star_ang_mom;
	
	double min_mass, max_mass;
	double gas_min_mass, gas_max_mass, dark_min_mass, dark_max_mass, star_min_mass, star_max_mass;
	double min_radius, max_radius;
	double gas_min_radius, gas_max_radius, dark_min_radius, dark_max_radius, star_min_radius, star_max_radius;
	double min_velocity, max_velocity;
	double gas_min_velocity, gas_max_velocity, dark_min_velocity, dark_max_velocity, star_min_velocity, star_max_velocity;
	double dark_min_eps, dark_max_eps, star_min_eps, star_max_eps;
	double min_phi, max_phi;
	double gas_min_phi, gas_max_phi, dark_min_phi, dark_max_phi, star_min_phi, star_max_phi;
	double gas_min_rho, gas_max_rho;
	double gas_min_temp, gas_max_temp;
	double gas_min_hsmooth, gas_max_hsmooth;
	double gas_min_metals, gas_max_metals, star_min_metals, star_max_metals;
	double star_min_tform, star_max_tform;
	
	void relocate_center_of_mass(const Vector& new_com);
	void set_center_of_mass_velocity(const Vector& new_com_vel);
	
	void outputStats(std::ostream& os);
};

vector<Real> readTipsyArray(std::istream& is);
vector<Vector3D<Real> > readTipsyVector(std::istream& is);

#endif //TIPSYFILE_H
