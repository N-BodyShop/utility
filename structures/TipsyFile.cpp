/** @file TipsyFile.cpp
 @author Graeme Lufkin (gwl@u.washington.edu)
 @date Created April 30, 2002
 @version 2.2
 */

#include <fstream>
#include <algorithm>

#include "TipsyFile.h"

namespace Tipsy {

using std::sort;
using std::unique;

//create a new tipsy file with the specified numbers of particles
TipsyFile::TipsyFile(const std::string& fn, int nGas, int nDark, int nStar) : native(true), success(false), filename(fn), h(nGas, nDark, nStar), gas(nGas), darks(nDark), stars(nStar) {
	success = true;
}

//open a tipsyfile from disk
TipsyFile::TipsyFile(const std::string& fn) : native(true), success(false), myReader(fn), filename(fn) {
	loadfile();
}

//open a tipsyfile from a stream
TipsyFile::TipsyFile(std::istream& is) : native(true), success(false), myReader(is), filename("") {
	loadfile();
}

bool TipsyFile::reload(const std::string& fn) {
	filename = fn;
	TipsyReader r(filename);
	myReader.takeOverStream(r);
	return loadfile();
}

bool TipsyFile::reload(std::istream& is) {
	filename = "";
	TipsyReader r(is);
	myReader.takeOverStream(r);
	return loadfile();
}

bool TipsyFile::saveAll() const {
	std::ofstream outfile(filename.c_str()); //open file
	if(!outfile)
		return false;
	bool result = saveAll(outfile);
	outfile.close();
	return result;
}

template <typename Iterator>
void write_vector(std::ostream& os, Iterator begin, Iterator end) {
	std::size_t itemSize = (*begin).sizeBytes;
	for(Iterator iter = begin; iter != end; ++iter) {
		os.write(reinterpret_cast<const char *>(&(*iter)), itemSize);
		if(!os)
			return;
	}
}

bool TipsyFile::saveAll(std::ostream& os) const {
	if(!os)
		return false;
	os.write(reinterpret_cast<const char *>(&h), header::sizeBytes); //write out the header
	if(!os)
		return false;
	
	//endian-ness check
	unsigned int bob = 3;
	unsigned char* c = reinterpret_cast<unsigned char *>(&bob);
	if(c[3] == bob) { //we're big-endian, write the pad
		bob = 0;
		os.write(reinterpret_cast<const char *>(&bob), 4);
		if(!os)
			return false;
	}
	
	write_vector(os, gas.begin(), gas.end());
	if(!os)
		return false;
	write_vector(os, darks.begin(), darks.end());
	if(!os)
		return false;
	write_vector(os, stars.begin(), stars.end());
	if(!os)
		return false;
	else
		return true;
}

bool TipsyFile::save() const {
	std::ofstream outfile(filename.c_str());
	if(!outfile)
		return false;
	bool result = save(outfile);
	outfile.close();
	return result;
}

bool TipsyFile::save(std::ostream& os) const {
	//sort the marked indices
	sort(markedGas.begin(), markedGas.end());
	sort(markedDarks.begin(), markedDarks.end());
	sort(markedStars.begin(), markedStars.end());
	//remove any duplicate entries
	markedGas.erase(unique(markedGas.begin(), markedGas.end()), markedGas.end());
	markedDarks.erase(unique(markedDarks.begin(), markedDarks.end()), markedDarks.end());
	markedStars.erase(unique(markedStars.begin(), markedStars.end()), markedStars.end());
	header newh = h;
	newh.nsph = gas.size() - markedGas.size();
	newh.ndark = darks.size() - markedDarks.size();
	newh.nstar = stars.size() - markedStars.size();
	newh.nbodies = newh.nsph + newh.ndark + newh.nstar;

	if(!os)
		return false;
	os.write(reinterpret_cast<const char *>(&newh), header::sizeBytes); //write out the header
	if(!os)
		return false;
	
	//endian-ness check
	unsigned int bob = 3;
	unsigned char* c = reinterpret_cast<unsigned char *>(&bob);
	if(c[3] == bob) { //we're big-endian, write the pad
		bob = 0;
		os.write(reinterpret_cast<const char *>(&bob), 4);
		if(!os)
			return false;
	}

	int n;
	if(markedGas.size() == 0) {
		write_vector(os, gas.begin(), gas.end());
		if(!os)
			return false;
	} else {
		std::vector<gas_particle>::const_iterator startGas = gas.begin();
		for(std::vector<int>::iterator iter = markedGas.begin(); iter != markedGas.end(); ++iter) {
			n = gas.begin() + *iter - startGas;
			write_vector(os, startGas, startGas + n);
			if(!os)
				return false;
			startGas += n + 1;
		}
		write_vector(os, startGas, gas.end());
	}

	if(markedDarks.size() == 0) {
		write_vector(os, darks.begin(), darks.end());
		if(!os)
			return false;
	} else {
		std::vector<dark_particle>::const_iterator startDark = darks.begin();
		write_vector(os, darks.begin(), darks.end());
		for(std::vector<int>::iterator iter = markedDarks.begin(); iter != markedDarks.end(); ++iter) {
			n = darks.begin() + *iter - startDark;
			write_vector(os, startDark, startDark + n);
			if(!os)
				return false;
			startDark += n + 1;
		}
		write_vector(os, startDark, darks.end());
	}

	if(markedStars.size() == 0) {
		write_vector(os, stars.begin(), stars.end());
		if(!os)
			return false;
	} else {
		std::vector<star_particle>::const_iterator startStar = stars.begin();
		for(std::vector<int>::iterator iter = markedStars.begin(); iter != markedStars.end(); ++iter) {
			n = stars.begin() + *iter - startStar;
			write_vector(os, startStar, startStar + n);
			if(!os)
				return false;
			startStar += n + 1;
		}
		write_vector(os, startStar, stars.end());
	}
	
	return true;
}

bool TipsyFile::writeIOrdFile() {
	return writeIOrdFile(filename + std::string(".iOrd"));
}

bool TipsyFile::writeIOrdFile(const std::string& fn) {
	std::ofstream outfile(fn.c_str());
	if(!outfile)
		return false;
	bool result = writeIOrdFile(outfile);
	outfile.close();
	return result;
}

bool TipsyFile::writeIOrdFile(std::ostream& os) {
	if(!os)
		return false;
	
	//sort the marked indices
	sort(markedGas.begin(), markedGas.end());
	sort(markedDarks.begin(), markedDarks.end());
	sort(markedStars.begin(), markedStars.end());
	//remove any duplicate entries
	markedGas.erase(unique(markedGas.begin(), markedGas.end()), markedGas.end());
	markedDarks.erase(unique(markedDarks.begin(), markedDarks.end()), markedDarks.end());
	markedStars.erase(unique(markedStars.begin(), markedStars.end()), markedStars.end());

	remakeHeader();
	const char* id = "iOrd";
	//write out an identifier
	os.write(id, 4 * sizeof(char));
	if(!os)
		return false;
	//write out the next available gas index
	os.write(reinterpret_cast<const char *>(&h.nsph), sizeof(int));
	if(!os)
		return false;
	//write out the next available dark index
	os.write(reinterpret_cast<const char *>(&h.ndark), sizeof(int));
	if(!os)
		return false;
	//write out the next available star index
	os.write(reinterpret_cast<const char *>(&h.nstar), sizeof(int));
	if(!os)
		return false;
	
	std::vector<int>::iterator markedIter;
	markedIter = markedGas.begin();
	for(int i = 0; i < h.nsph; i++) {
		if(markedIter == markedGas.end()) {
			os.write(reinterpret_cast<const char *>(&i), sizeof(int));
			if(!os)
				return false;
		} else if(i == *markedIter)
			++markedIter;
		else {
			os.write(reinterpret_cast<const char *>(&i), sizeof(int));
			if(!os)
				return false;
		}
	}
	markedIter = markedDarks.begin();
	for(int i = 0; i < h.ndark; i++) {
		if(markedIter == markedDarks.end()) {
			os.write(reinterpret_cast<const char *>(&i), sizeof(int));
			if(!os)
				return false;
		} else if(i == *markedIter)
			++markedIter;
		else {
			os.write(reinterpret_cast<const char *>(&i), sizeof(int));
			if(!os)
				return false;
		}
	}
	markedIter = markedStars.begin();
	for(int i = 0; i < h.nstar; i++) {
		if(markedIter == markedStars.end()) {
			os.write(reinterpret_cast<const char *>(&i), sizeof(int));
			if(!os)
				return false;
		} else if(i == *markedIter)
			++markedIter;
		else {
			os.write(reinterpret_cast<const char *>(&i), sizeof(int));
			if(!os)
				return false;
		}
	}
	return true;
}

bool TipsyFile::loadfile() {
	success = false;
	
	if(!myReader.status())
		return false;
	h = myReader.getHeader();
	native = myReader.isNative();
	
	gas.resize(h.nsph);
	for(int i = 0; i < h.nsph; ++i) {
		if(!myReader.getNextGasParticle(gas[i]))
			return false;
		boundingBox.grow(gas[i].pos);
	}
	darks.resize(h.ndark);
	for(int i = 0; i < h.ndark; ++i) {
		if(!myReader.getNextDarkParticle(darks[i]))
			return false;
		boundingBox.grow(darks[i].pos);
	}
	stars.resize(h.nstar);
	for(int i = 0; i < h.nstar; ++i) {
		if(!myReader.getNextStarParticle(stars[i]))
			return false;
		boundingBox.grow(stars[i].pos);
	}
	
	success = true;
	return true;
}

//read part of a tipsy file from a file
PartialTipsyFile::PartialTipsyFile(const std::string& fn, int begin, int end) : myReader(fn), filename(fn), beginParticle(begin), endParticle(end) {
	loadPartial();
}

//read part of a tipsy file from a stream
PartialTipsyFile::PartialTipsyFile(std::istream& is, int begin, int end) : myReader(is), filename(""), beginParticle(begin), endParticle(end) {
	loadPartial();
}

bool PartialTipsyFile::reload(const std::string& fn) {
	filename = fn;
	TipsyReader r(filename);
	myReader.takeOverStream(r);
	return loadPartial();
}

bool PartialTipsyFile::reload(std::istream& is) {
	filename = "";
	TipsyReader r(is);
	myReader.takeOverStream(r);
	return loadPartial();
}

bool PartialTipsyFile::loadPartial() {
	success = false;
	
	if(!myReader.status())
		return false;
	
	fullHeader = myReader.getHeader();
	native = myReader.isNative();
	
	//check that the range of particles asked for is reasonable
	if(beginParticle < 0 || beginParticle > fullHeader.nbodies) {
		return false;
	}
	
	//just read all the rest of the particles if endParticle is bogus
	if(endParticle <= beginParticle || endParticle > fullHeader.nbodies)
		endParticle = fullHeader.nbodies;
	
	//we will need to seek at least this far in the file
	//int seekPosition = headerSizeBytes + (native ? 0 : sizeof(int));
	
	//fill in the header for our part
	h.nbodies = endParticle - beginParticle;
	h.ndim = fullHeader.ndim;
	h.time = fullHeader.time;
	//figure out the left and right boundaries of each type of particle
	h.nsph = std::max(std::min(fullHeader.nsph, endParticle) - beginParticle, 0);
	h.ndark = std::max(std::min(endParticle, fullHeader.nsph + fullHeader.ndark) - std::max(fullHeader.nsph, beginParticle), 0);
	h.nstar = std::max(endParticle - std::max(fullHeader.nsph + fullHeader.ndark, beginParticle), 0);
	
	if(!myReader.seekParticleNum(beginParticle))
		return false;
	
	gas.resize(h.nsph);
	for(int i = 0; i < h.nsph; ++i) {
		if(!myReader.getNextGasParticle(gas[i]))
			return false;
	}
	darks.resize(h.ndark);
	for(int i = 0; i < h.ndark; ++i) {
		if(!myReader.getNextDarkParticle(darks[i]))
			return false;
	}
	stars.resize(h.nstar);
	for(int i = 0; i < h.nstar; ++i) {
		if(!myReader.getNextStarParticle(stars[i]))
			return false;
	}

	success = true;
	return true;
}

TipsyStats::TipsyStats(TipsyFile* tfile) {
	tf = tfile;
	
	total_mass = gas_mass = dark_mass = star_mass = 0;
	
	volume = density = 0;
	
	//set initial value for the bounding box with the first particle we can find
	if(tf->h.nsph > 0)
		bounding_box = OrientedBox<double>(Vector(tf->gas[0].pos), Vector(tf->gas[0].pos));
	else if(tf->h.ndark > 0)
		bounding_box = OrientedBox<double>(tf->darks[0].pos, tf->darks[0].pos);
	else if(tf->h.nstar > 0)
		bounding_box = OrientedBox<double>(tf->stars[0].pos, tf->stars[0].pos);
	else //no particles in the file, so why did you ask for stats?!
		return;
	
	center = size = Vector(0, 0, 0);
	
	center_of_mass = gas_com = dark_com = star_com = Vector(0, 0, 0);
	center_of_mass_velocity = gas_com_vel = dark_com_vel = star_com_vel = Vector(0, 0, 0);
	angular_momentum = gas_ang_mom = dark_ang_mom = star_ang_mom = Vector(0, 0, 0);

	min_mass = gas_min_mass = dark_min_mass = star_min_mass = min_radius \
			= gas_min_radius = dark_min_radius = star_min_radius \
			= min_velocity = gas_min_velocity = dark_min_velocity \
			= star_min_velocity = dark_min_eps = star_min_eps = min_phi \
			= gas_min_phi = dark_min_phi = star_min_phi = gas_min_rho \
			= gas_min_temp = gas_min_hsmooth  = gas_min_metals = star_min_metals \
			= star_min_tform = HUGE_VAL;
	max_mass = gas_max_mass = dark_max_mass = star_max_mass = max_radius \
			= gas_max_radius = dark_max_radius = star_max_radius \
			= max_velocity = gas_max_velocity = dark_max_velocity \
			= star_max_velocity = dark_max_eps = star_max_eps = max_phi \
			= gas_max_phi = dark_max_phi = star_max_phi = gas_max_rho \
			= gas_max_temp = gas_max_hsmooth  = gas_max_metals = star_max_metals \
			= star_max_tform = -HUGE_VAL;
	
	int i;
	Vector position, velocity;

	//collect stats on gas particles
	for(i = 0; i < tf->h.nsph; ++i) {
		gas_mass += tf->gas[i].mass;
		position = tf->gas[i].pos;
		velocity = tf->gas[i].vel;
		bounding_box.grow(position);
		
		gas_com += tf->gas[i].mass * position;
		gas_com_vel += tf->gas[i].mass * velocity;
		gas_ang_mom += tf->gas[i].mass * cross(position, velocity);
		
		if(tf->gas[i].mass < gas_min_mass)
			gas_min_mass = tf->gas[i].mass;
		if(tf->gas[i].mass > gas_max_mass)
			gas_max_mass = tf->gas[i].mass;
		if(position.length() < gas_min_radius)
			gas_min_radius = position.length();
		if(position.length() > gas_max_radius)
			gas_max_radius = position.length();
		if(velocity.length() < gas_min_velocity)
			gas_min_velocity = velocity.length();
		if(velocity.length() > gas_max_velocity)
			gas_max_velocity = velocity.length();
		if(tf->gas[i].phi < gas_min_phi)
			gas_min_phi = tf->gas[i].phi;
		if(tf->gas[i].phi > gas_max_phi)
			gas_max_phi = tf->gas[i].phi;
		if(tf->gas[i].hsmooth < gas_min_hsmooth)
			gas_min_hsmooth = tf->gas[i].hsmooth;
		if(tf->gas[i].hsmooth > gas_max_hsmooth)
			gas_max_hsmooth = tf->gas[i].hsmooth;
		if(tf->gas[i].rho < gas_min_rho)
			gas_min_rho = tf->gas[i].rho;
		if(tf->gas[i].rho > gas_max_rho)
			gas_max_rho = tf->gas[i].rho;
		if(tf->gas[i].temp < gas_min_temp)
			gas_min_temp = tf->gas[i].temp;
		if(tf->gas[i].temp > gas_max_temp)
			gas_max_temp = tf->gas[i].temp;
		if(tf->gas[i].metals < gas_min_metals)
			gas_min_metals = tf->gas[i].metals;
		if(tf->gas[i].metals > gas_max_metals)
			gas_max_metals = tf->gas[i].metals;
	}
	
	//collect stats on dark particles
	for(i = 0; i < tf->h.ndark; ++i) {
		dark_mass += tf->darks[i].mass;
		position = tf->darks[i].pos;
		velocity = tf->darks[i].vel;
		bounding_box.grow(position);
		
		dark_com += tf->darks[i].mass * position;
		dark_com_vel += tf->darks[i].mass * velocity;
		dark_ang_mom += tf->darks[i].mass * cross(position, velocity);
		
		if(tf->darks[i].mass < dark_min_mass)
			dark_min_mass = tf->darks[i].mass;
		if(tf->darks[i].mass > dark_max_mass)
			dark_max_mass = tf->darks[i].mass;
		if(position.length() < dark_min_radius)
			dark_min_radius = position.length();
		if(position.length() > dark_max_radius)
			dark_max_radius = position.length();
		if(velocity.length() < dark_min_velocity)
			dark_min_velocity = velocity.length();
		if(velocity.length() > dark_max_velocity)
			dark_max_velocity = velocity.length();
		if(tf->darks[i].phi < dark_min_phi)
			dark_min_phi = tf->darks[i].phi;
		if(tf->darks[i].phi > dark_max_phi)
			dark_max_phi = tf->darks[i].phi;
		if(tf->darks[i].eps < dark_min_eps)
			dark_min_eps = tf->darks[i].eps;
		if(tf->darks[i].eps > dark_max_eps)
			dark_max_eps = tf->darks[i].eps;
	}

	//collect stats on star particles
	for(i = 0; i < tf->h.nstar; ++i) {
		star_mass += tf->stars[i].mass;
		position = tf->stars[i].pos;
		velocity = tf->stars[i].vel;
		bounding_box.grow(position);
		
		star_com += tf->stars[i].mass * position;
		star_com_vel += tf->stars[i].mass * velocity;
		star_ang_mom += tf->stars[i].mass * cross(position, velocity);
		
		if(tf->stars[i].mass < star_min_mass)
			star_min_mass = tf->stars[i].mass;
		if(tf->stars[i].mass > star_max_mass)
			star_max_mass = tf->stars[i].mass;
		if(position.length() < star_min_radius)
			star_min_radius = position.length();
		if(position.length() > star_max_radius)
			star_max_radius = position.length();
		if(velocity.length() < star_min_velocity)
			star_min_velocity = velocity.length();
		if(velocity.length() > star_max_velocity)
			star_max_velocity = velocity.length();
		if(tf->stars[i].phi < star_min_phi)
			star_min_phi = tf->stars[i].phi;
		if(tf->stars[i].phi > star_max_phi)
			star_max_phi = tf->stars[i].phi;
		if(tf->stars[i].eps < star_min_eps)
			star_min_eps = tf->stars[i].eps;
		if(tf->stars[i].eps > star_max_eps)
			star_max_eps = tf->stars[i].eps;
		if(tf->stars[i].metals < star_min_metals)
			star_min_metals = tf->stars[i].metals;
		if(tf->stars[i].metals > star_max_metals)
			star_max_metals = tf->stars[i].metals;
		if(tf->stars[i].tform < star_min_tform)
			star_min_tform = tf->stars[i].tform;
		if(tf->stars[i].tform > star_max_tform)
			star_max_tform = tf->stars[i].tform;
	}
	
	//get the total mass of the system
	total_mass = gas_mass + dark_mass + star_mass;
	
	//get the total center of mass position and velocity
	if(total_mass != 0) {
		center_of_mass = (gas_com + dark_com + star_com) / total_mass;
		center_of_mass_velocity = (gas_com_vel + dark_com_vel + star_com_vel) / total_mass;
	}
	
	//get total angular momentum about the origin
	angular_momentum = gas_ang_mom + dark_ang_mom + star_ang_mom;
	//make specific (divide by total mass)
	angular_momentum /= total_mass;
	//subtract off the angular momentum of the center of mass
	angular_momentum -= cross(center_of_mass, center_of_mass_velocity);

	//sort out the individual centers of mass and angular momenta
	if(gas_mass != 0) {
		gas_com /= gas_mass;
		gas_com_vel /= gas_mass;
		gas_ang_mom /= gas_mass;
		gas_ang_mom -= cross(gas_com, gas_com_vel);
	}
	if(dark_mass != 0) {
		dark_com /= dark_mass;
		dark_com_vel /= dark_mass;
		dark_ang_mom /= dark_mass;
		dark_ang_mom -= cross(dark_com, dark_com_vel);
	}
	if(star_mass != 0) {
		star_com /= star_mass;
		star_com_vel /= star_mass;
		star_ang_mom /= star_mass;
		star_ang_mom -= cross(star_com, star_com_vel);
	}
		
	//get position stats
	center = bounding_box.center();
	size = bounding_box.greater_corner - bounding_box.lesser_corner;
	volume = bounding_box.volume();
	if(volume > 0)
		density = total_mass / volume;
	
	//figure out aggregate ranges
	min_mass = (gas_min_mass < dark_min_mass ? gas_min_mass : dark_min_mass);
	min_mass = (min_mass < star_min_mass ? min_mass : star_min_mass);
	max_mass = (gas_max_mass > dark_max_mass ? gas_max_mass : dark_max_mass);
	max_mass = (max_mass > star_max_mass ? max_mass : star_max_mass);
	min_radius = (gas_min_radius < dark_min_radius ? gas_min_radius : dark_min_radius);
	min_radius = (min_radius < star_min_radius ? min_radius : star_min_radius);
	max_radius = (gas_max_radius > dark_max_radius ? gas_max_radius : dark_max_radius);
	max_radius = (max_radius > star_max_radius ? max_radius : star_max_radius);
	min_velocity = (gas_min_velocity < dark_min_velocity ? gas_min_velocity : dark_min_velocity);
	min_velocity = (min_velocity < star_min_velocity ? min_velocity : star_min_velocity);
	max_velocity = (gas_max_velocity > dark_max_velocity ? gas_max_velocity : dark_max_velocity);
	max_velocity = (max_velocity > star_max_velocity ? max_velocity : star_max_velocity);
	min_phi = (gas_min_phi < dark_min_phi ? gas_min_phi : dark_min_phi);
	min_phi = (min_phi < star_min_phi ? min_phi : star_min_phi);
	max_phi = (gas_max_phi > dark_max_phi ? gas_max_phi : dark_max_phi);
	max_phi = (max_phi > star_max_phi ? max_phi : star_max_phi);	
}

void TipsyStats::outputStats(std::ostream& os) {
	using std::endl;
	os << "Info on tipsy file " << tf->filename << endl;
	if(tf->isNative())
		os << "x86 file format (little-endian)" << endl;
	else
		os << "SGI file format (big-endian)" << endl;
	os << "Simulation time: " << tf->h.time << endl;
	os << "NBodies: " << tf->h.nbodies << " (" << tf->h.nsph << " gas, " << tf->h.ndark << " dark, " << tf->h.nstar << " star)" << endl;
	os << "Bounding box: " << bounding_box << endl;
	os << "Box center: " << center << endl;
	os << "Box volume: " << volume << " Box size: " << size << endl;
	os << "Density: " << density << " Total mass: " << total_mass << endl;
	os << "Total center of mass: " << center_of_mass << endl;
	os << "Total center of mass velocity: " << center_of_mass_velocity << endl;
	os << "Total angular momentum: " << angular_momentum << endl;
	
	if(tf->h.nsph) {
		os << "Gas Stats: (" << tf->h.nsph << " particles)" << endl;
		os << "Gas mass: " << gas_mass << endl;
		os << "Gas center of mass: " << gas_com << endl;
		os << "Gas center of mass velocity: " << gas_com_vel << endl;
		os << "Gas angular momentum: " << gas_ang_mom << endl;
		os << "Range of gas masses: (" << gas_min_mass << " <=> " << gas_max_mass << ")" << endl;
		os << "Range of gas radii: (" << gas_min_radius << " <=> " << gas_max_radius << ")" << endl;
		os << "Range of gas velocities: (" << gas_min_velocity << " <=> " << gas_max_velocity << ")" << endl;
		os << "Range of gas rho: (" << gas_min_rho << " <=> " << gas_max_rho << ")" << endl;
		os << "Range of gas temperature: (" << gas_min_temp << " <=> " << gas_max_temp << ")" << endl;
		os << "Range of gas hsmooth: (" << gas_min_hsmooth << " <=> " << gas_max_hsmooth << ")" << endl;
		os << "Range of gas metals: (" << gas_min_metals << " <=> " << gas_max_metals << ")" << endl;
		os << "Range of gas phi: (" << gas_min_phi << " <=> " << gas_max_phi << ")" << endl;
	} else
		os << "No gas particles." << endl;
	
	if(tf->h.ndark) {
		os << "Dark Stats: (" << tf->h.ndark << " particles)" << endl;
		os << "Dark mass: " << dark_mass << endl;
		os << "Dark center of mass: " << dark_com << endl;
		os << "Dark center of mass velocity: " << dark_com_vel << endl;
		os << "Dark angular momentum: " << dark_ang_mom << endl;
		os << "Range of dark masses: (" << dark_min_mass << " <=> " << dark_max_mass << ")" << endl;
		os << "Range of dark radii: (" << dark_min_radius << " <=> " << dark_max_radius << ")" << endl;
		os << "Range of dark velocities: (" << dark_min_velocity << " <=> " << dark_max_velocity << ")" << endl;
		os << "Range of eps: (" << dark_min_eps << " <=> " << dark_max_eps << ")" << endl;
		os << "Range of dark phi: (" << dark_min_phi << " <=> " << dark_max_phi << endl;
	} else
		os << "No dark particles." << endl;
	
	if(tf->h.nstar) {
		os << "Star Stats: (" << tf->h.nstar << " particles)" << endl;
		os << "Star mass: " << star_mass << endl;
		os << "Star center of mass: " << star_com << endl;
		os << "Star center of mass velocity: " << star_com_vel << endl;
		os << "Star angular momentum: " << star_ang_mom << endl;
		os << "Range of star masses: (" << star_min_mass << " <=> " << star_max_mass << ")" << endl;
		os << "Range of star radii: (" << star_min_radius << " <=> " << star_max_radius << ")" << endl;
		os << "Range of star velocities: (" << star_min_velocity << " <=> " << star_max_velocity << ")" << endl;
		os << "Range of star metals: (" << star_min_metals << " <=> " << star_max_metals << ")" << endl;
		os << "Range of star formation time: (" << star_min_tform << " <=> " << star_max_tform << ")" << endl;
		os << "Range of star eps: (" << star_min_eps << " <=> " << star_max_eps << ")" << endl;
		os << "Range of star phi: (" << star_min_phi << " <=> " << star_max_phi << ")" << endl;
	} else
		os << "No star particles." << endl;
	
	os << "Aggregate Stats:" << endl;
	os << "Range of masses: (" << min_mass << " <=> " << max_mass << ")" << endl;
	os << "Range of radii: (" << min_radius << " <=> " << max_radius << ")" << endl;
	os << "Range of velocities: (" << min_velocity << " <=> " << max_velocity << ")" << endl;
	os << "Range of phi: (" << min_phi << " <=> " << max_phi << ")" << endl;
}

void TipsyStats::relocate_center_of_mass(const Vector& new_com) {
	int i;
	for(i = 0; i < tf->h.nsph; ++i)
		tf->gas[i].pos -= new_com;
	
	for(i = 0; i < tf->h.ndark; ++i)
		tf->darks[i].pos -= new_com;
	
	for(i = 0; i < tf->h.nstar; ++i)
		tf->stars[i].pos -= new_com;
}

void TipsyStats::set_center_of_mass_velocity(const Vector& new_com_vel) {
	int i;
	for(i = 0; i < tf->h.nsph; ++i)
		tf->gas[i].vel -= new_com_vel;
	
	for(i = 0; i < tf->h.ndark; ++i)
		tf->darks[i].vel -= new_com_vel;

	for(i = 0; i < tf->h.nstar; ++i)
		tf->stars[i].vel -= new_com_vel;
}

std::vector<Real> readTipsyArray(std::istream& is) {
	int num;
	is >> num;
	std::vector<Real> bad;
	if(!is || (num <= 0))
		return bad;
	std::vector<Real> arrayvals(num);
	for(int i = 0; i < num; ++i) {
		is >> arrayvals[i];
		if(!is)
			return bad;
	}
	return arrayvals;
}

std::vector<Vector3D<Real> > readTipsyVector(std::istream& is) {
	int num;
	is >> num;
	std::vector<Vector3D<Real> > bad;
	if(!is || (num <= 0))
		return bad;
	std::vector<Vector3D<Real> > vectorvals(num);
	int i;
	for(i = 0; i < num; ++i) {
		is >> vectorvals[i].x;
		if(!is)
			return bad;
	}
	for(i = 0; i < num; ++i) {
		is >> vectorvals[i].y;
		if(!is)
			return bad;
	}
	for(i = 0; i < num; ++i) {
		is >> vectorvals[i].z;
		if(!is)
			return bad;
	}
	return vectorvals;
}

} //close namespace Tipsy
