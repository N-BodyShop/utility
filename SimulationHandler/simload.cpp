//simload.cpp

#include <iostream>
#include <string>
#include <list>
#include <vector>

#include "Simulation.h"
#include "SiXFormat.h"

using namespace SimulationHandling;
using namespace std;

extern "C" void CkRegisterMainModule(void) {

}

int main(int argc, char** argv) {

	if(argc < 2) {
		cerr << "Usage: simload simulation_directory" << endl;
		return 1;
	}
	
	Simulation* sim = new SiXFormatReader(argv[1]);
	cout << "Simulation " << sim->name << " contains " << sim->size() << " families" << endl;
	
	for(Simulation::iterator iter = sim->begin(); iter != sim->end(); ++iter) {
		cout << "Family " << iter->first << " has the following attributes:" << endl;
		ParticleFamily& family = iter->second;
		for(AttributeMap::iterator attrIter = family.attributes.begin(); attrIter != family.attributes.end(); ++attrIter)
			cout << attrIter->first << " scalar min: " << getScalarMin(attrIter->second) << " max: " << getScalarMax(attrIter->second) << endl;
	}
	
	sim->loadAttribute("dark", "position");
	sim->loadAttribute("dark", "mass");
	
	ParticleFamily& darks = sim->getFamily("dark");
	Vector3D<float>* positions = darks.getAttribute<Vector3D<float> >("position");
	float* masses = darks.getAttribute<float>("mass");
	if(!masses)
		cerr << "Didn't get masses!" << endl;
	if(!positions)
		cerr << "Didn't get positions!" << endl;
	if(masses && positions) {
		Vector3D<double> cm;
		double totalMass = 0;
		u_int64_t N = darks.count.totalNumParticles;
		cout << "Finding cm for " << N << " dark particles" << endl;
		for(u_int64_t i = 0; i < N; ++i) {
			cm += masses[i] * positions[i];
			totalMass += masses[i];
		}
		cm /= totalMass;
		cout << "Total dark mass: " << totalMass << endl;
		cout << "Dark center of mass: " << cm << endl;
	} else
		cerr << "Couldn't load attributes successfully!" << endl;
	
	sim->release();
	
	delete sim;
	
	cerr << "Done." << endl;
}
