/** @file ParticleGroup.h
 @author Graeme Lufkin (gwl@u.washington.edu)
 @date Created November 21, 2003
 @version 1.0
 */
 
#ifndef PARTICLEGROUP_H__ujiomre9yq5389879861rwqeuhmi
#define PARTICLEGROUP_H__ujiomre9yq5389879861rwqeuhmi

#include <vector>

#include "Simulation.h"

namespace SimulationHandling {

class ParticleGroup : public std::vector<u_int64_t> {
	bool status;
	
public:
	
	//typedef std::vector<u_int64_t>::iterator iterator;

	template <typename T>
	ParticleGroup(const ParticleFamily& family, const std::string& attributeName, T minValue, T maxValue) : status(false) {
		if(T* values = family.getAttribute(attributeName, Type2Type<T>())) {
			for(u_int64_t i = 0; i < family.count.numParticles; ++i) {
				if(minValue <= values[i] && values[i] <= maxValue)
					push_back(i);
			}
		}
	}
	
	template <typename T>
	ParticleGroup(const ParticleFamily& family, const std::string& attributeName, Vector3D<T> minValue, Vector3D<T> maxValue) : status(false) {
		if(Vector3D<T>* values = family.getAttribute(attributeName, Type2Type<Vector3D<T> >())) {
			OrientedBox<T> box(minValue, maxValue);
			for(u_int64_t i = 0; i < family.count.numParticles; ++i) {
				if(box.contains(values[i]))
					push_back(i);
			}
		}
	}
	
	
	bool getStatus() const {
		return status;
	}
};


template <typename T>
bool createGroup(const std::string& groupName, const std::string& attributeName, T minVal, T maxVal) {
	
}

} //close namespace SimulationHandling

#endif //PARTICLEGROUP_H__ujiomre9yq5389879861rwqeuhmi
