/** \file TipsyFile.cpp
 \author Graeme Lufkin (gwl@u.washington.edu)
 \date Created April 30, 2002
 \version 2.0
 \todo Fully document this file
 */

#include "TipsyFile.h"

namespace Tipsy {

//create a new tipsy file with the specified numbers of particles
TipsyFile::TipsyFile(const string& fn, int nGas, int nDark = 0, int nStar = 0) : native(true), success(false), filename(fn), h(nGas, nDark, nStar), gas(nGas), darks(nDark), stars(nStar) {
	success = true;
}

//open a tipsyfile from disk
TipsyFile::TipsyFile(const string& fn) : native(true), success(false), filename(fn) {
	reload(filename);
}

//open a tipsyfile from a stream
TipsyFile::TipsyFile(std::istream& is) : native(true), success(false), filename("") {
	reload(is);
}

void TipsyFile::reload(const std::string& fn) {
	filename = fn;
	ifstream infile(filename.c_str());
	loadfile(infile);
	infile.close();
}

void TipsyFile::reload(std::istream& is) {
	loadfile(is);
}

//process xdr versions of the data structures used in tipsy files
int TipsyFile::xdr_convert_header(header &h) {
	return (xdr_double(&xdrs, &h.time) 
		&& xdr_int(&xdrs, &h.nbodies)
		&& xdr_int(&xdrs, &h.ndim) 
		&& xdr_int(&xdrs, &h.nsph)
		&& xdr_int(&xdrs, &h.ndark)
		&& xdr_int(&xdrs, &h.nstar));
}

int TipsyFile::xdr_convert_gas(gas_particle &g) {
	return (xdr_float(&xdrs, &g.mass)
		&& xdr_float(&xdrs, &g.pos[0])
		&& xdr_float(&xdrs, &g.pos[1])
		&& xdr_float(&xdrs, &g.pos[2])
		&& xdr_float(&xdrs, &g.vel[0])
		&& xdr_float(&xdrs, &g.vel[1])
		&& xdr_float(&xdrs, &g.vel[2])
		&& xdr_float(&xdrs, &g.rho)
		&& xdr_float(&xdrs, &g.temp)
		&& xdr_float(&xdrs, &g.hsmooth)
		&& xdr_float(&xdrs, &g.metals)
		&& xdr_float(&xdrs, &g.phi));
}

int TipsyFile::xdr_convert_dark(dark_particle &d) {
	return (xdr_float(&xdrs, &d.mass)
		&& xdr_float(&xdrs, &d.pos[0])
		&& xdr_float(&xdrs, &d.pos[1])
		&& xdr_float(&xdrs, &d.pos[2])
		&& xdr_float(&xdrs, &d.vel[0])
		&& xdr_float(&xdrs, &d.vel[1])
		&& xdr_float(&xdrs, &d.vel[2])
		&& xdr_float(&xdrs, &d.eps)
		&& xdr_float(&xdrs, &d.phi));
}

int TipsyFile::xdr_convert_star(star_particle &s) {
	return (xdr_float(&xdrs, &s.mass)
		&& xdr_float(&xdrs, &s.pos[0])
		&& xdr_float(&xdrs, &s.pos[1])
		&& xdr_float(&xdrs, &s.pos[2])
		&& xdr_float(&xdrs, &s.vel[0])
		&& xdr_float(&xdrs, &s.vel[1])
		&& xdr_float(&xdrs, &s.vel[2])
		&& xdr_float(&xdrs, &s.metals)
		&& xdr_float(&xdrs, &s.tform)
		&& xdr_float(&xdrs, &s.eps)
		&& xdr_float(&xdrs, &s.phi));
}

//write a TipsyFile to file
void TipsyFile::save() const {
	std::ofstream outfile(filename.c_str()); //open file
	save(outfile);
	outfile.close();
}

//write a tipsy file to stream
void TipsyFile::save(std::ostream& os) const {
	os.write(reinterpret_cast<const char *>(&h), sizeof(header)); //write out the header
	os.write(reinterpret_cast<const char *>(gas.begin()), h.nsph * sizeof(gas_particle)); //write out the particles
	os.write(reinterpret_cast<const char *>(darks.begin()), h.ndark * sizeof(dark_particle));
	os.write(reinterpret_cast<const char *>(stars.begin()), h.nstar * sizeof(star_particle));
}

void TipsyFile::loadfile(std::istream& in) {
	success = false;
	
	//read in the header
	in.read(&h, sizeof(header));
	
	if(h.ndim != MAXDIM) { //check for validity
		//try xdr format
		xdrmem_create(&xdrs, (char *) &h, sizeof(header), XDR_DECODE);
		if(!xdr_convert_header(h) || h.ndim != MAXDIM) { //wasn't xdr format either			
			cerr << "Tried to load file " << filename << " which had crazy dimension.\ntipsy file not loaded." << endl;
			h.nbodies = h.nsph = h.ndark = h.nstar = 0;
			return;
		}
		native = false;
		//xdr format has an integer pad in the header, which we don't need, but must skip
		int pad = 0;
		in.read((char *) &pad, sizeof(int));
	}

	//check and allocate memory for gas particles
	if(h.nsph > 0) {
		gas.resize(h.nsph);
		in.read(gas.begin(), h.nsph * sizeof(gas_particle)); //read in the gas particles
		if(!native) { //do xdr processing
			//create an xdr stream in memory
			xdrmem_create(&xdrs, reinterpret_cast<char *>(gas.begin()), h.nsph * sizeof(gas_particle), XDR_DECODE);
			for(vector<gas_particle>::iterator iter = gas.begin(), end = gas.end(); iter != end; ++iter)
				if(!xdr_convert_gas(*iter)) //convert each particle and check for validity
					return;
		}
	}

	//check and allocate memory for dark matter particles
	if(h.ndark > 0) {
		darks.resize(h.ndark);
		in.read(darks.begin(), h.ndark * sizeof(dark_particle));
		if(!native) {
			xdrmem_create(&xdrs, reinterpret_cast<char *>(darks.begin()), h.ndark * sizeof(dark_particle), XDR_DECODE);
			for(vector<dark_particle>::iterator iter = darks.begin(), end = darks.end(); iter != end; ++iter)
				if(!xdr_convert_dark(*iter)) //convert each particle and check for validity
					return;
		}
	}

	//check and allocate memory for star particles
	if(h.nstar > 0) {
		stars.resize(h.nstar);
		in.read(stars.begin(), h.nstar * sizeof(star_particle));
		if(!native) {
			xdrmem_create(&xdrs, reinterpret_cast<char *>(stars.begin()), h.nstar * sizeof(star_particle), XDR_DECODE);
			for(vector<star_particle>::iterator iter = stars.begin(), end = stars.end(); iter != end; ++iter)
				if(!xdr_convert_star(*iter)) //convert each particle and check for validity
					return;
		}
	}

	success = true;
}

//read part of a tipsy file from a file
PartialTipsyFile::PartialTipsyFile(const string& fn, int begin = 0, int end = 1) : beginParticle(begin), endParticle(end) {
	reload(fn);
}

//read part of a tipsy file from a stream
PartialTipsyFile::PartialTipsyFile(std::istream& is, int begin = 0, int end = 1) : beginParticle(begin), endParticle(end) {
	reload(is);
}

void PartialTipsyFile::reload(const std::string& fn) {
	filename = fn;
	std::ifstream infile(filename.c_str());
	if(!infile)
		return;
	loadPartial(infile);
	infile.close();	
}

void PartialTipsyFile::reload(std::istream& is) {
	loadPartial(is);
}

void PartialTipsyFile::loadPartial(std::istream& in) {
	success = false;
	
	//read in the header
	in.read(&fullHeader, sizeof(header));
	
	if(fullHeader.ndim != MAXDIM) { //check for validity
		//try xdr format
		xdrmem_create(&xdrs, reinterpret_cast<char *>(&fullHeader), sizeof(header), XDR_DECODE);
		if(!xdr_convert_header(fullHeader) || fullHeader.ndim != MAXDIM) { //wasn't xdr format either			
			cerr << "Tried to load file " << filename << " which had crazy dimension.\ntipsy file not loaded." << endl;
			h.nbodies = h.nsph = h.ndark = h.nstar = 0;
			fullHeader.nbodies = fullHeader.nsph = fullHeader.ndark = fullHeader.nstar = 0;
			return;
		}
		native = false;
		//xdr format has an integer pad in the header, which we don't need, but must skip
		int pad = 0;
		in.read((char *) &pad, sizeof(int));
	}
	
	if(beginParticle < 0 || beginParticle > fullHeader.nbodies) {
		cerr << "You asked for an unreasonable portion of the tipsy file" << endl;
		return;
	}
	
	//just read all the rest of the particles if endParticle is bogus
	if(endParticle <= beginParticle || endParticle > fullHeader.nbodies)
		endParticle = fullHeader.nbodies;
	
	//we will need to seek at least this far in the file
	int seekPosition = sizeof(header);
	
	h.nbodies = endParticle - beginParticle;
	h.ndim = fullHeader.ndim;
	h.time = fullHeader.time;
	//figure out the left and right boundaries of each type of particle
	h.nsph = max(min(fullHeader.nsph, endParticle) - beginParticle, 0);
	h.ndark = max(min(endParticle, fullHeader.nsph + fullHeader.ndark) - max(fullHeader.nsph, beginParticle), 0);
	h.nstar = max(endParticle - max(fullHeader.nsph + fullHeader.ndark, beginParticle), 0);

	//determine how many bytes we need to seek
	if(beginParticle < fullHeader.nsph)
		seekPosition += beginParticle * sizeof(gas_particle);
	else if(beginParticle < fullHeader.nsph + fullHeader.ndark)
		seekPosition += fullHeader.nsph * sizeof(gas_particle) + (beginParticle - fullHeader.nsph) * sizeof(dark_particle);
	else
		seekPosition += fullHeader.nsph * sizeof(gas_particle) + fullHeader.ndark * sizeof(dark_particle) + (beginParticle - fullHeader.nsph - fullHeader.ndark) * sizeof(star_particle);
	
	in.seekg(seekPosition);
	
	//check and allocate memory for gas
	if(h.nsph > 0) {
		gas.resize(h.nsph);
		in.read(gas.begin(), h.nsph * sizeof(gas_particle)); //read in the gas particles
		if(!native) { //do xdr processing
			//create an xdr stream in memory
			xdrmem_create(&xdrs, reinterpret_cast<char *>(gas.begin()), h.nsph * sizeof(gas_particle), XDR_DECODE);
			for(vector<gas_particle>::iterator iter = gas.begin(), end = gas.end(); iter != end; ++iter)
				if(!xdr_convert_gas(*iter)) //convert each particle and check for validity
					return;
		}
	}

	//check and allocate memory for dark matter
	if(h.ndark > 0) {
		darks.resize(h.ndark);
		in.read(darks.begin(), h.ndark * sizeof(dark_particle));
		if(!native) {
			xdrmem_create(&xdrs, reinterpret_cast<char *>(darks.begin()), h.ndark * sizeof(dark_particle), XDR_DECODE);
			for(vector<dark_particle>::iterator iter = darks.begin(), end = darks.end(); iter != end; ++iter)
				if(!xdr_convert_dark(*iter)) //convert each particle and check for validity
					return;
		}
	}

	//check and allocate memory for stars
	if(h.nstar > 0) {
		stars.resize(h.nstar);
		in.read(stars.begin(), h.nstar * sizeof(star_particle));
		if(!native) {
			xdrmem_create(&xdrs, reinterpret_cast<char *>(stars.begin()), h.nstar * sizeof(star_particle), XDR_DECODE);
			for(vector<star_particle>::iterator iter = stars.begin(), end = stars.end(); iter != end; ++iter)
				if(!xdr_convert_star(*iter)) //convert each particle and check for validity
					return;
		}
	}

	success = true;
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
			= star_min_tform = 1E1024;
	max_mass = gas_max_mass = dark_max_mass = star_max_mass = max_radius \
			= gas_max_radius = dark_max_radius = star_max_radius \
			= max_velocity = gas_max_velocity = dark_max_velocity \
			= star_max_velocity = dark_max_eps = star_max_eps = max_phi \
			= gas_max_phi = dark_max_phi = star_max_phi = gas_max_rho \
			= gas_max_temp = gas_max_hsmooth  = gas_max_metals = star_max_metals \
			= star_max_tform = -1E1024;
	
	int i;
	Vector position, velocity;

	//collect stats on gas particles
	for(i = 0; i < tf->h.nsph; i++) {
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
	for(i = 0; i < tf->h.ndark; i++) {
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
	for(i = 0; i < tf->h.nstar; i++) {
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
	for(i = 0; i < tf->h.nsph; i++)
		tf->gas[i].pos -= new_com;
	
	for(i = 0; i < tf->h.ndark; i++)
		tf->darks[i].pos -= new_com;
	
	for(i = 0; i < tf->h.nstar; i++)
		tf->stars[i].pos -= new_com;
}

void TipsyStats::set_center_of_mass_velocity(const Vector& new_com_vel) {
	int i;
	for(i = 0; i < tf->h.nsph; i++)
		tf->gas[i].vel -= new_com_vel;
	
	for(i = 0; i < tf->h.ndark; i++)
		tf->darks[i].vel -= new_com_vel;

	for(i = 0; i < tf->h.nstar; i++)
		tf->stars[i].vel -= new_com_vel;
}

vector<Real> readTipsyArray(std::istream& is) {
	int num;
	is >> num;
	vector<Real> bad;
	if(!is || (num <= 0))
		return bad;
	vector<Real> arrayvals(num);
	for(int i = 0; i < num; i++) {
		is >> arrayvals[i];
		if(!is)
			return bad;
	}
	return arrayvals;
}

vector<Vector3D<Real> > readTipsyVector(std::istream& is) {
	int num;
	is >> num;
	vector<Vector3D<Real> > bad;
	if(!is || (num <= 0))
		return bad;
	vector<Vector3D<Real> > vectorvals(num);
	int i;
	for(i = 0; i < num; i++) {
		is >> vectorvals[i].x;
		if(!is)
			return bad;
	}
	for(i = 0; i < num; i++) {
		is >> vectorvals[i].y;
		if(!is)
			return bad;
	}
	for(i = 0; i < num; i++) {
		is >> vectorvals[i].z;
		if(!is)
			return bad;
	}
	return vectorvals;
}

} //close namespace Tipsy
