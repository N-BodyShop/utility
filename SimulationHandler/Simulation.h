/** @file Simulation.h
 @author Graeme Lufkin (gwl@u.washington.edu)
 @date Created September 17, 2003
 @version 1.0
 */
 
#ifndef SIMULATION_H
#define SIMULATION_H

#include <string>
#include <map>
#include <list>

#include <sys/types.h>

#include "tree_xdr.h"

namespace SimulationHandling {

using namespace TypeHandling;

struct ParticleCount {
	u_int64_t numParticles;
	u_int64_t startParticle;
	u_int64_t totalNumParticles;
	
	ParticleCount() : numParticles(0), startParticle(0), totalNumParticles(0) { }
};

typedef std::map<std::string, TypedArray> AttributeMap;

class ParticleFamily {
public:

	std::string familyName;
	ParticleCount count;
	AttributeMap attributes;
	
	ParticleFamily(const std::string& name = "") : familyName(name) { }
	
	template <typename T>
	T* getAttribute(const std::string& attributeName) {
		AttributeMap::iterator iter = attributes.find(attributeName);
		if(iter == attributes.end())
			return 0;
		return iter->second.getArray(Type2Type<T>());
	}
	
	void addAttribute(const std::string& attributeName, const TypeHandling::TypedArray& arr) {
		if(arr.length == count.numParticles)
			attributes[attributeName] = arr;
	}
	
	template <typename T>
	void addAttribute(const std::string& attributeName, T* data) {
		TypeHandling::TypedArray& arr = attributes[attributeName];
		arr.dimensions = Type2Dimensions<T>::dimensions;
		arr.code = Type2Code<T>::code;
		arr.length = count.numParticles;
		arr.data = data;
		arr.calculateMinMax();
	}
	
	bool releaseAttribute(const std::string& attributeName) {
		AttributeMap::iterator iter = attributes.find(attributeName);
		if(iter == attributes.end())
			return false;
		iter->second.release();
		return true;
	}
	
	void release() {
		for(AttributeMap::iterator iter = attributes.begin(); iter != attributes.end(); ++iter)
			iter->second.release();
	}
};

class Simulation : public std::map<std::string, ParticleFamily> {
public:
		
	std::string name;
	
	ParticleFamily& getFamily(const std::string& familyName) {
		return operator[](familyName);
	}
	
	/// A sub-class implements the loading of values for a particular attribute
	virtual bool loadAttribute(const std::string& familyName, const std::string& attributeName, u_int64_t numParticles = 0, const u_int64_t startParticle = 0) = 0;

	void release() {
		for(iterator iter = begin(); iter != end(); ++iter)
			iter->second.release();
	}
	
	u_int64_t numParticles() const {
		u_int64_t numParticles = 0;
		for(const_iterator iter = begin(); iter != end(); ++iter)
			numParticles += iter->second.count.numParticles;
		return numParticles;
	}
	
	u_int64_t totalNumParticles() const {
		u_int64_t numParticles = 0;
		for(const_iterator iter = begin(); iter != end(); ++iter)
			numParticles += iter->second.count.totalNumParticles;
		return numParticles;
	}
	
	virtual ~Simulation() { }
	
};

class SimulationWriter {
public:
	/** Write the given simulation to disk in a format determined by the implementing sub-class.
	All non-null attribute arrays will be written to disk.  To prevent an attribute
	from being written do disk, explicitly release it before calling this method. */
	virtual bool save(const Simulation*, const std::string&) = 0;

	//virtual bool write(const Simulation*, const std::string&) = 0;
};

} //close namespace SimulationHandling

#endif //SIMULATION_H
