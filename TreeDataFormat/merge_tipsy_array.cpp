/** @file merge_tipsy_array.cpp
 Given a SiX format file and a Tipsy array file, merge the array
 values into the SiX file.
 */

#include <iostream>
#include <fstream>

#include "Simulation.h"
#include "SiXFormat.h"

template <typename IndexType, typename LengthType, typename ValueType>
void reorder_array(IndexType* indexArray, LengthType N, ValueType* valueArray) {
	IndexType index;
	for(LengthType i = 0; i < N; ++i) {
		index = i;
		while((index = indexArray[index]) < i);
		swap(valueArray[i], valueArray[index]);
	}
}

using namespace std;
using namespace SimulationHandling;

bool mergeAttribute(Simulation* sim, const string& familyName, const string& attributeName, u_int64_t& numPrevious, istream& infile) {
	Simulation::iterator iter = sim->find(familyName);
	if(iter != sim->end()) { //there are some gas particles
		ParticleFamily& family = iter->second;
		sim->loadAttribute(iter->first, "uid", family.count.totalNumParticles);
		unsigned int* uids = family.getAttribute<unsigned int>("uid");
		if(!uids) {
			cerr << "SiX format did not contain Tipsy-order uids!" << endl;
			return false;
		}
		float* values = new float[family.count.totalNumParticles];
		for(unsigned int i = 0; i < family.count.totalNumParticles; ++i) {
			infile >> values[i];
			if(!infile) {
				cerr << "Problem reading values from array file" << endl;
				return false;
			}
		}
		
		reorder_array(uids, family.count.totalNumParticles, values);
		numPrevious += family.count.totalNumParticles;	
		family.releaseAttribute("uid");
		family.addAttribute(attributeName, values);
	}
	
	return true;
}

int main(int argc, char** argv) {
	
	if(argc < 4) {
		cerr << "Usage: merge_tipsy_array SiXFile arrayfile attribute_name" << endl;
		return 1;
	}
	
	Simulation* sim = new SiXFormatReader(argv[1]);
	if(sim->size() == 0) {
		cerr << "Problems opening SiX format file" << endl;
		return 2;
	}
	
	ifstream infile(argv[2]);
	if(!infile) {
		cerr << "Problems reading Tipsy array file" << endl;
		return 2;
	}
	
	int numParticles;
	//read number of particles from array file
	infile >> numParticles;
	if(!infile) {
		cerr << "Couldn't read number of particles from array file" << endl;
		return 2;
	}
	
	if(numParticles != sim->totalNumParticles()) {
		cerr << "Number of particles in array doesn't match SiX format file" << endl;
		return 2;
	}
	
	string attributeName = argv[3];
	
	u_int64_t numPrevious = 0;
	
	//merge the attribute into the Tipsy families
	if(!mergeAttribute(sim, "gas", attributeName, numPrevious, infile))
		return 3;
	if(!mergeAttribute(sim, "dark", attributeName, numPrevious, infile))
		return 3;
	if(!mergeAttribute(sim, "star", attributeName, numPrevious, infile))
		return 3;
	
	cerr << "Attribute merged from array file, writing to disk ..." << endl;
	
	SimulationWriter* writer = new SiXFormatWriter;
	writer->save(sim, "");
	
	delete writer;
	delete sim;
	
	cerr << "Done." << endl;
}
