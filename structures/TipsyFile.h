/** @file TipsyFile.h
 @author Graeme Lufkin (gwl@u.washington.edu)
 @date Created April 30, 2002
 @version 2.2
 */
 
#ifndef TIPSYFILE_H
#define TIPSYFILE_H

#include <iostream>
#include <vector>
#include <string>

#include "TipsyParticles.h"
#include "Vector3D.h"
#include "OrientedBox.h"

namespace Tipsy {

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
	
	header(int nGas = 0, int nDark = 0, int nStar = 0) : time(0), nbodies(nGas + nDark + nStar), ndim(MAXDIM), nsph(nGas), ndark(nDark), nstar(nStar) { }
	
	/// Output operator, used for formatted display
	friend std::ostream& operator<< (std::ostream& os, const header& h) {
		return os << "Time: " << h.time
			<< "\nnBodies: " << h.nbodies
			<< "\nnDim: " << h.ndim
			<< "\nnSPH: " << h.nsph
			<< "\nnDark: " << h.ndark
			<< "\nnStar: " << h.nstar;
	}
};

/** This class represents a tipsy format file in memory. */
class TipsyFile {
private:

	/// Load the file from disk
	bool loadfile(std::istream& in);
	
	mutable std::vector<int> markedGas;
	mutable std::vector<int> markedDarks;
	mutable std::vector<int> markedStars;
	
	bool native;

	bool success;
	
public:
	/// The filename of this tipsy file
	std::string filename;

	/// The header of this tipsy file
	header h;

	/// The array of gas particles
	std::vector<gas_particle> gas;
	/// The array of dark matter particles
	std::vector<dark_particle> darks;
	/// The array of star particles
	std::vector<star_particle> stars;
	
	/// The axis-aligned box that contains all the particles
	OrientedBox<Real> boundingBox;
	
	TipsyFile() : native(true), success(false) { }

	/// Create a blank tipsy file with specified number of particles
	TipsyFile(const std::string& fn, int nGas, int nDark = 0, int nStar = 0);
	
	/// Load from a file
	TipsyFile(const std::string& fn);
	
	/// Load from a stream
	TipsyFile(std::istream& is);

	/// Copy constructor
	TipsyFile(const TipsyFile& tf) : markedGas(tf.markedGas), markedDarks(tf.markedDarks), markedStars(tf.markedStars), native(tf.native), success(tf.success), filename(tf.filename), h(tf.h), gas(tf.gas), darks(tf.darks), stars(tf.stars), boundingBox(tf.boundingBox) { }
	
	~TipsyFile() { }
	
	TipsyFile& operator=(const TipsyFile& tf) {
		markedGas = tf.markedGas;
		markedDarks = tf.markedDarks;
		markedStars = tf.markedStars;
		native = tf.native;
		success = tf.success;
		filename = tf.filename;
		h = tf.h;
		gas = tf.gas;
		darks = tf.darks;
		stars = tf.stars;
		boundingBox = tf.boundingBox;
		return *this;
	}
	
	/// Load a new file into this object, destroying the original one
	bool reload(const std::string& fn);
	bool reload(std::istream& is);
	
	/// Save the current state to the given filename, not writing out particles marked for deletion
	bool save() const;
	/// Save the current state to a stream
	bool save(std::ostream& os) const;
	
	bool writeIOrdFile();
	bool writeIOrdFile(const std::string& fn);
	bool writeIOrdFile(std::ostream& os);
	
	/// Save the current state, including particles marked for deletion
	bool saveAll() const;
	bool saveAll(std::ostream& os) const;
	
	void addGasParticle(const gas_particle& p = gas_particle()) {
		gas.push_back(p);
		boundingBox.grow(p.pos);
		remakeHeader();
	}
	
	void addDarkParticle(const dark_particle& p = dark_particle()) {
		darks.push_back(p);
		boundingBox.grow(p.pos);
		remakeHeader();
	}
	
	void addStarParticle(const star_particle& p = star_particle()) {
		stars.push_back(p);
		boundingBox.grow(p.pos);
		remakeHeader();
	}
	
	bool markGasForDeletion(const int i) {
		if(i < 0 || i >= h.nsph)
			return false;
		markedGas.push_back(i);
		return true;
	}
	
	bool markDarkForDeletion(const int i) {
		if(i < 0 || i >= h.ndark)
			return false;
		markedDarks.push_back(i);
		return true;
	}

	bool markStarForDeletion(const int i) {
		if(i < 0 || i >= h.nstar)
			return false;
		markedStars.push_back(i);
		return true;
	}
	
	/// Is this file stored in native byte-order?
	bool isNative() const { return native; }

	/// Did the file load successfully?
	bool loadedSuccessfully() const { return success; }
	
	bool operator!() const { return !success; }
	
	void remakeHeader() {
		h.nsph = gas.size();
		h.ndark = darks.size();
		h.nstar = stars.size();
		h.nbodies = h.nsph + h.ndark + h.nstar;
	}
};

/** This class represents part of a tipsy file loaded into memory as if it were
 the whole thing.  It can be used in cases of low memory, or to split up one file
 among several processes.  It is essentially a read-only structure.  i.e., you
 cannot save, add particles, or remove particles from a PartialTipsyFile.  */
class PartialTipsyFile {
private:
		
	bool loadPartial(std::istream& is);
	bool native;
	bool success;
	
public:
	
	/// The filename of this tipsy file
	std::string filename;

	/// The header for the full file
	header fullHeader;
	/// The header for the part of the file we hold
	header h;

	/// The array of gas particles
	std::vector<gas_particle> gas;
	/// The array of dark matter particles
	std::vector<dark_particle> darks;
	/// The array of star particles
	std::vector<star_particle> stars;
	
	/// The axis-aligned box that contains all the particles
	OrientedBox<Real> boundingBox;
	
	/// The on-disk index of the first particle we hold
	int beginParticle;
	/// The on-disk index of one past the last particle we hold
	int endParticle;
	
	PartialTipsyFile() : native(true), success(false), beginParticle(0), endParticle(0) { }
	PartialTipsyFile(const std::string& fn, int begin = 0, int end = 1);
	PartialTipsyFile(std::istream& is, int begin = 0, int end = 1);
	PartialTipsyFile(const PartialTipsyFile& ptf) : native(ptf.native), success(ptf.success), filename(ptf.filename), fullHeader(ptf.fullHeader), h(ptf.h), gas(ptf.gas), darks(ptf.darks), stars(ptf.stars), beginParticle(ptf.beginParticle), endParticle(ptf.endParticle) { }
		
	PartialTipsyFile& operator=(const PartialTipsyFile& ptf) {
		native = ptf.native;
		success = ptf.success;
		filename = ptf.filename;
		fullHeader = ptf.fullHeader;
		h = ptf.h;
		gas = ptf.gas;
		darks = ptf.darks;
		stars = ptf.stars;
		beginParticle = ptf.beginParticle;
		endParticle = ptf.endParticle;
		return *this;
	}
	
	//reloading a partial file
	bool reload(const std::string& fn);
	bool reload(std::istream& is);
	
	/// Is this file stored in native byte-order?
	bool isNative() const { return native; }

	/// Did the file load successfully?
	bool loadedSuccessfully() const { return success; }
	
	bool operator!() const { return !success; }
};

/** This class hold statistics calculated from a tipsy file. */
class TipsyStats {
private:
	TipsyFile* tf;
public:
	TipsyStats() { }
	TipsyStats(TipsyFile* tfile);
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

std::vector<Real> readTipsyArray(std::istream& is);
std::vector<Vector3D<Real> > readTipsyVector(std::istream& is);

} //close namespace Tipsy

#endif //TIPSYFILE_H
