/** @file Group.h
 @author Graeme Lufkin (gwl@u.washington.edu)
 @date Created June 28, 2004
 @version 1.0
 */
 
#ifndef GROUP_H__7h237h327yh4th7y378w37aw7ho4wvta8
#define GROUP_H__7h237h327yh4th7y378w37aw7ho4wvta8

#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/iterator/iterator_facade.hpp>

#include "Simulation.h"

namespace SimulationHandling {

using namespace TypeHandling;

/** The base class for group iterators.
 Pointers to instances of this class or subclasses form the implementation
 of the concrete \c GroupIterator class.
 Individual types of group will subclass this and implement their
 own increment, equal and dereference methods.
 */
class GroupIteratorInstance {
	u_int64_t index;
public:
	
	GroupIteratorInstance() : index(0) { }
	GroupIteratorInstance(u_int64_t start) : index(start) { }
	virtual ~GroupIteratorInstance() { }
	
	/// For criterion-based groups, this method should move to the next member of the parent group that satisfies the criterion, or the end
	virtual void increment() {
		++index;
	}
	
	/// For non-trivial groups, check both array index and criterion data
	virtual bool equal(GroupIteratorInstance* const& other) const {
		return index == other->index;
	}
	
	/// This should return the array index of the current group member
	virtual u_int64_t dereference() const {
		return index;
	}
};

/** A concrete class that uses the pimpl idiom to handle many types of
 group.  Instances of \c Group and its subclasses will provide
 objects of this type holding the appropriate \c impl pointer
 to the actual iterator object.
 */
class GroupIterator : public boost::iterator_facade<GroupIterator, u_int64_t, boost::forward_traversal_tag, u_int64_t> {
public:
	
	GroupIterator() { }

	explicit GroupIterator(boost::shared_ptr<GroupIteratorInstance> impl_) : impl(impl_) { }
	
	void setImplementation(boost::shared_ptr<GroupIteratorInstance> impl_) {
		impl = impl_;
	}
	
private:
		
	friend class boost::iterator_core_access;

	void increment() {
		if(impl)
			impl->increment();
	}
	
	bool equal(GroupIterator const& other) const {
		return (impl.get() ? (other.impl.get() ? impl->equal(other.impl.get()) : false) : false);
	}
	
	u_int64_t dereference() const {
		return (impl.get() ? impl->dereference() : 0);
	}

	boost::shared_ptr<GroupIteratorInstance> impl;
};

/** A \c Group object has virtual functions that provide iterators for
 the particular kind of group, and lists the families it is valid for.
 The big picture is: keep a collection of pointers to \c Group objects.
 When processing a group, use the virtual member functions to create
 concrete \c GroupIterator objects.  These iterators use the pimpl
 idiom to implement the appropriate behavior for the group they work for.
 */
class Group {
	//disable copying
	Group(Group const&);
	Group& operator=(Group const&);
	
public:
	
	Group(boost::shared_ptr<Group> const& parent = boost::shared_ptr<Group>()) : parentGroup(parent) { }
	virtual ~Group() { }
		
	typedef std::vector<std::string> GroupFamilies;
	GroupFamilies families;
	
	virtual GroupIterator make_begin_iterator(std::string const& familyName) = 0;
	virtual GroupIterator make_end_iterator(std::string const& familyName) = 0;

protected:
	boost::shared_ptr<Group> parentGroup;
};

/// The "All" group iterator is the default behavior, provide a typedef.
typedef GroupIteratorInstance AllGroupIterator;

/// The All group doesn't have a parent
class AllGroup : public Group {
	Simulation const& sim;
public:
		
	AllGroup(Simulation const& s) : sim(s) {
		families.reserve(sim.size());
		for(Simulation::const_iterator iter = sim.begin(); iter != sim.end(); ++iter)
			families.push_back(iter->first);
	}
	
	GroupIterator make_begin_iterator(std::string const& familyName) {
		boost::shared_ptr<AllGroupIterator> p(new AllGroupIterator(0));
		return GroupIterator(p);
	}

	GroupIterator make_end_iterator(std::string const& familyName) {
		Simulation::const_iterator simIter = sim.find(familyName);
		if(simIter == sim.end())
			return make_begin_iterator(familyName);
		boost::shared_ptr<AllGroupIterator> p(new AllGroupIterator(simIter->second.count.numParticles));
		return GroupIterator(p);
	}
};

template <typename T>
class AttributeRangeIterator : public GroupIteratorInstance {
	template <typename>
	friend class AttributeRangeGroup;
	
	//u_int64_t index;
	//parent group idiom: Instead of index, use iterator into parent group
	GroupIterator parentIter;
	T minValue;
	T maxValue;
	T const* array; //bare pointer is okay here, default copy constructor will do what we want
	u_int64_t length;
	
public:
	
	AttributeRangeIterator(GroupIterator parentBegin, T minValue_, T maxValue_, T const* array_, u_int64_t length_) : parentIter(parentBegin), minValue(minValue_), maxValue(maxValue_), array(array_), length(length_) {
		while(*parentIter < length && (array[*parentIter] < minValue || maxValue < array[*parentIter]))
			++parentIter;
	}

	void increment() {
		//parent group idiom: increment parentIter instead of index
		if(*parentIter < length)
			for(++parentIter; *parentIter < length && (array[*parentIter] < minValue || maxValue < array[*parentIter]); ++parentIter);
	}
	
	bool equal(GroupIteratorInstance* const& other) const {
		if(AttributeRangeIterator* const o = dynamic_cast<AttributeRangeIterator* const>(other))
			return minValue == o->minValue && maxValue == o->maxValue && array == o->array && parentIter == o->parentIter;
		else
			return false;
	}
	
	u_int64_t dereference() const {
		//parent group idiom: instead of index value, dereference iterator
		return *parentIter;
	}
};

template <typename T>
class AttributeRangeGroup : public Group {
	Simulation const& sim;
	std::string attributeName;
	T minValue;
	T maxValue;
public:
		
	AttributeRangeGroup(Simulation const& s, boost::shared_ptr<Group> const& parent, std::string const& attributeName_, T minValue_, T maxValue_) : Group(parent), sim(s), attributeName(attributeName_), minValue(minValue_), maxValue(maxValue_) {
		for(Simulation::const_iterator simIter = sim.begin(); simIter != sim.end(); ++simIter) {
			AttributeMap::const_iterator attrIter = simIter->second.attributes.find(attributeName);
			if(attrIter != simIter->second.attributes.end())
				families.push_back(simIter->first);
		}
	}
	
	GroupIterator make_begin_iterator(std::string const& familyName) {
		Simulation::const_iterator simIter = sim.find(familyName);
		if(simIter == sim.end())
			return make_end_iterator(familyName);
		TypedArray const& array = simIter->second.attributes.find(attributeName)->second;
		//parent group idiom: start with parent group's begin iterator
		GroupIterator parentBegin = parentGroup->make_begin_iterator(familyName);
		boost::shared_ptr<AttributeRangeIterator<T> > p(new AttributeRangeIterator<T>(parentBegin, minValue, maxValue, array.getArray(Type2Type<T>()), array.length));
		return GroupIterator(p);
	}

	GroupIterator make_end_iterator(std::string const& familyName) {
		//parent group idiom: parent group can handle end iterator
		return parentGroup->make_end_iterator(familyName);
	}
};

boost::shared_ptr<Group> make_AttributeRangeGroup(Simulation const& sim, boost::shared_ptr<Group> const& parent, std::string const& attributeName, double minValue, double maxValue) {
	boost::shared_ptr<Group> p;
	for(Simulation::const_iterator simIter = sim.begin(); simIter != sim.end(); ++simIter) {
		AttributeMap::const_iterator attrIter = simIter->second.attributes.find(attributeName);
		if(attrIter != simIter->second.attributes.end()) {
			//found a family with the attribute
			TypedArray const& arr = attrIter->second;
			if(arr.dimensions == 1) {
				switch(arr.code) {
					case int8:
						p.reset(new AttributeRangeGroup<Code2Type<int8>::type>(sim, parent, attributeName, static_cast<Code2Type<int8>::type>(minValue), static_cast<Code2Type<int8>::type>(maxValue)));
						break;
					case uint8:
						p.reset(new AttributeRangeGroup<Code2Type<uint8>::type>(sim, parent, attributeName, static_cast<Code2Type<uint8>::type>(minValue), static_cast<Code2Type<uint8>::type>(maxValue)));
						break;
					case int16:
						p.reset(new AttributeRangeGroup<Code2Type<int16>::type>(sim, parent, attributeName, static_cast<Code2Type<int16>::type>(minValue), static_cast<Code2Type<int16>::type>(maxValue)));
						break;
					case uint16:
						p.reset(new AttributeRangeGroup<Code2Type<uint16>::type>(sim, parent, attributeName, static_cast<Code2Type<uint16>::type>(minValue), static_cast<Code2Type<uint16>::type>(maxValue)));
						break;
					case int32:
						p.reset(new AttributeRangeGroup<Code2Type<int32>::type>(sim, parent, attributeName, static_cast<Code2Type<int32>::type>(minValue), static_cast<Code2Type<int32>::type>(maxValue)));
						break;
					case uint32:
						p.reset(new AttributeRangeGroup<Code2Type<uint32>::type>(sim, parent, attributeName, static_cast<Code2Type<uint32>::type>(minValue), static_cast<Code2Type<uint32>::type>(maxValue)));
						break;
					case int64:
						p.reset(new AttributeRangeGroup<Code2Type<int64>::type>(sim, parent, attributeName, static_cast<Code2Type<int64>::type>(minValue), static_cast<Code2Type<int64>::type>(maxValue)));
						break;
					case uint64:
						p.reset(new AttributeRangeGroup<Code2Type<uint64>::type>(sim, parent, attributeName, static_cast<Code2Type<uint64>::type>(minValue), static_cast<Code2Type<uint64>::type>(maxValue)));
						break;
					case float32:
						p.reset(new AttributeRangeGroup<Code2Type<float32>::type>(sim, parent, attributeName, static_cast<Code2Type<float32>::type>(minValue), static_cast<Code2Type<float32>::type>(maxValue)));
						break;
					case float64:
						p.reset(new AttributeRangeGroup<Code2Type<float64>::type>(sim, parent, attributeName, static_cast<Code2Type<float64>::type>(minValue), static_cast<Code2Type<float64>::type>(maxValue)));
						break;
				}
			}
			break;
		}
	}
	return p;
}

/*
class NoCriterionGroupIterator : public GroupIteratorInstance {
	GroupIterator parentIter;	
public:
	
	NoCriterionGroupIterator(GroupIterator parentBegin) : parentIter(parentBegin) { }

	void increment() {
		if(*parentIter < length)
			++parentIter;
	}
	
	bool equal(GroupIteratorInstance* const& other) const {
		if(NoCriterionGroupIterator* const o = dynamic_cast<NoCriterionGroupIterator* const>(other))
			return parentIter == o->parentIter;
		else
			return false;
	}
	
	u_int64_t dereference() const {
		return *parentIter;
	}
};
*/
		
template <typename T>
class SphericalIterator : public GroupIteratorInstance {
	template <typename>
	friend class SphericalGroup;
	
	u_int64_t index;
	Vector3D<T> centerVector;
	T radiusValue;
	Vector3D<T> const* array; //bare pointer is okay here, default copy constructor will do what we want
	u_int64_t length;
	
	void resetBegin() {
		for(index = 0; index < length && ((array[index]-centerVector).length() > radiusValue); ++index);
	}

public:
	
	SphericalIterator(u_int64_t start, Vector3D<T> centerVector_, T radiusValue_, Vector3D<T> const* array_, u_int64_t length_) : index(start), centerVector(centerVector_), radiusValue(radiusValue_), array(array_), length(length_) { }

	void increment() {
		if(index >= length - 1)
			index = length;
		else
			for(++index; index < length && ((array[index]-centerVector).length() > radiusValue); ++index);
	}
	
	bool equal(GroupIteratorInstance* const& other) const {
		if(SphericalIterator* const o = dynamic_cast<SphericalIterator* const>(other))
			return centerVector == o->centerVector && radiusValue == o->radiusValue && array == o->array && index == o->index;
		else
			return false;
	}
	
	u_int64_t dereference() const {
		return index;
	}
};

template <typename T>
class SphericalGroup : public Group {
	Simulation const& sim;
	std::string attributeName;
	Vector3D<T> centerVector;
	T radiusValue;
public:
		
	SphericalGroup(Simulation const& s, std::string const& attributeName_, Vector3D<T> centerVector_, T radiusValue_) : sim(s), attributeName(attributeName_), centerVector(centerVector_), radiusValue(radiusValue_) {
		for(Simulation::const_iterator simIter = sim.begin(); simIter != sim.end(); ++simIter) {
			AttributeMap::const_iterator attrIter = simIter->second.attributes.find(attributeName);
			if(attrIter != simIter->second.attributes.end())
				families.push_back(simIter->first);
		}
	}
	
	GroupIterator make_begin_iterator(std::string const& familyName) {
		boost::shared_ptr<SphericalIterator<T> > p;
		Simulation::const_iterator simIter = sim.find(familyName);
		if(simIter == sim.end()) {
		    Vector3D<T> tmp(0,0,0);
		    p.reset(new SphericalIterator<T>(0, tmp, 0, 0, 0));
		    }
		else {
			AttributeMap::const_iterator attrIter = simIter->second.attributes.find(attributeName);
			TypedArray const& array = attrIter->second;
			p.reset(new SphericalIterator<T>(0, centerVector, radiusValue, array.getArray(Type2Type<Vector3D<T> >()), array.length));
			p->resetBegin();
		}
		return GroupIterator(p);
	}

	GroupIterator make_end_iterator(std::string const& familyName) {
		boost::shared_ptr<SphericalIterator<T> > p;
		Simulation::const_iterator simIter = sim.find(familyName);
		if(simIter == sim.end()) {
		    Vector3D<T> tmp(0,0,0);
		    p.reset(new SphericalIterator<T>(0, tmp, 0, 0, 0));
		    }
		else {
			AttributeMap::const_iterator attrIter = simIter->second.attributes.find(attributeName);
			TypedArray const& array = attrIter->second;
			p.reset(new SphericalIterator<T>(array.length, centerVector, radiusValue, array.getArray(Type2Type<Vector3D<T> >()), array.length));
		}
		return GroupIterator(p);
	}
};

static
boost::shared_ptr<Group> make_SphericalGroup(Simulation const& sim, std::string const& attributeName, Vector3D<double> centerVector, double radiusValue) {
	boost::shared_ptr<Group> p;
	for(Simulation::const_iterator simIter = sim.begin(); simIter != sim.end(); ++simIter) {
		AttributeMap::const_iterator attrIter = simIter->second.attributes.find(attributeName);
		if(attrIter != simIter->second.attributes.end()) {
			//found a family with the attribute
			TypedArray const& arr = attrIter->second;
			if(arr.dimensions == 1) {
				switch(arr.code) {
					case int8:
						p.reset(new SphericalGroup<Code2Type<int8>::type>(sim, attributeName, static_cast<Vector3D<Code2Type<int8>::type> >(centerVector), static_cast<Code2Type<int8>::type>(radiusValue)));
						break;
					case uint8:
						p.reset(new SphericalGroup<Code2Type<uint8>::type>(sim, attributeName, static_cast<Vector3D<Code2Type<uint8>::type> >(centerVector), static_cast<Code2Type<uint8>::type>(radiusValue)));
						break;
					case int16:
						p.reset(new SphericalGroup<Code2Type<int16>::type>(sim, attributeName, static_cast<Vector3D<Code2Type<int16>::type> >(centerVector), static_cast<Code2Type<int16>::type>(radiusValue)));
						break;
					case uint16:
						p.reset(new SphericalGroup<Code2Type<uint16>::type>(sim, attributeName, static_cast<Vector3D<Code2Type<uint16>::type> >(centerVector), static_cast<Code2Type<uint16>::type>(radiusValue)));
						break;
					case int32:
						p.reset(new SphericalGroup<Code2Type<int32>::type>(sim, attributeName, static_cast<Vector3D<Code2Type<int32>::type> >(centerVector), static_cast<Code2Type<int32>::type>(radiusValue)));
						break;
					case uint32:
						p.reset(new SphericalGroup<Code2Type<uint32>::type>(sim, attributeName, static_cast<Vector3D<Code2Type<uint32>::type> >(centerVector), static_cast<Code2Type<uint32>::type>(radiusValue)));
						break;
					case int64:
						p.reset(new SphericalGroup<Code2Type<int64>::type>(sim, attributeName, static_cast<Vector3D<Code2Type<int64>::type> >(centerVector), static_cast<Code2Type<int64>::type>(radiusValue)));
						break;
					case uint64:
						p.reset(new SphericalGroup<Code2Type<uint64>::type>(sim, attributeName, static_cast<Vector3D<Code2Type<uint64>::type> >(centerVector), static_cast<Code2Type<uint64>::type>(radiusValue)));
						break;
					case float32:
						p.reset(new SphericalGroup<Code2Type<float32>::type>(sim, attributeName, static_cast<Vector3D<Code2Type<float32>::type> >(centerVector), static_cast<Code2Type<float32>::type>(radiusValue)));
						break;
					case float64:
						p.reset(new SphericalGroup<Code2Type<float64>::type>(sim, attributeName, static_cast<Vector3D<Code2Type<float64>::type> >(centerVector), static_cast<Code2Type<float64>::type>(radiusValue)));
						break;
				}
			}
			break;
		}
	}
	return p;
}

class FamilyGroup : public Group {
	Simulation const& sim;
public:
		
	FamilyGroup(Simulation const& s, boost::shared_ptr<Group> const& parent, std::string const& familyName) : Group(parent), sim(s) {
		families.push_back(familyName);
	}
	
	GroupIterator make_begin_iterator(std::string const& familyName) {
		Simulation::const_iterator simIter = sim.find(familyName);
		if(simIter == sim.end())
			return make_end_iterator(familyName);
		return parentGroup->make_begin_iterator(familyName);
	}

	GroupIterator make_end_iterator(std::string const& familyName) {
		return parentGroup->make_end_iterator(familyName);
	}
};

//Static groups

//Tree-based groups

} //close namespace SimulationHandling

#endif //GROUP_H__7h237h327yh4th7y378w37aw7ho4wvta8
