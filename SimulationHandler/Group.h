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

	virtual void increment() {
		++index;
	}
	
	virtual bool equal(GroupIteratorInstance* const& other) const {
		return index == other->index;
	}
	
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
	
	Group() { }
	virtual ~Group() { }
		
	typedef std::vector<std::string> GroupFamilies;
	GroupFamilies families;
	
	virtual GroupIterator make_begin_iterator(std::string const& familyName) = 0;
	virtual GroupIterator make_end_iterator(std::string const& familyName) = 0;
		
};

/// The "All" group iterator is the default behavior, provide a typedef.
typedef GroupIteratorInstance AllGroupIterator;

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
	
	u_int64_t index;
	T minValue;
	T maxValue;
	T const* array; //bare pointer is okay here, default copy constructor will do what we want
	u_int64_t length;
	
	void resetBegin() {
		for(index = 0; index < length && (array[index] < minValue || maxValue < array[index]); ++index);
	}

public:
	
	AttributeRangeIterator(u_int64_t start, T minValue_, T maxValue_, T const* array_, u_int64_t length_) : index(start), minValue(minValue_), maxValue(maxValue_), array(array_), length(length_) { }

	void increment() {
		if(index >= length - 1)
			index = length;
		else
			for(++index; index < length && (array[index] < minValue || maxValue < array[index]); ++index);
	}
	
	bool equal(GroupIteratorInstance* const& other) const {
		if(AttributeRangeIterator* const o = dynamic_cast<AttributeRangeIterator* const>(other))
			return minValue == o->minValue && maxValue == o->maxValue && array == o->array && index == o->index;
		else
			return false;
	}
	
	u_int64_t dereference() const {
		return index;
	}
};

template <typename T>
class AttributeRangeGroup : public Group {
	Simulation const& sim;
	std::string attributeName;
	T minValue;
	T maxValue;
public:
		
	AttributeRangeGroup(Simulation const& s, std::string const& attributeName_, T minValue_, T maxValue_) : sim(s), attributeName(attributeName_), minValue(minValue_), maxValue(maxValue_) {
		for(Simulation::const_iterator simIter = sim.begin(); simIter != sim.end(); ++simIter) {
			AttributeMap::const_iterator attrIter = simIter->second.attributes.find(attributeName);
			if(attrIter != simIter->second.attributes.end())
				families.push_back(simIter->first);
		}
	}
	
	GroupIterator make_begin_iterator(std::string const& familyName) {
		boost::shared_ptr<AttributeRangeIterator<T> > p;
		Simulation::const_iterator simIter = sim.find(familyName);
		if(simIter == sim.end())
			p.reset(new AttributeRangeIterator<T>(0, 0, 0, 0, 0));
		else {
			AttributeMap::const_iterator attrIter = simIter->second.attributes.find(attributeName);
			TypedArray const& array = attrIter->second;
			p.reset(new AttributeRangeIterator<T>(0, minValue, maxValue, array.getArray(Type2Type<T>()), array.length));
			p->resetBegin();
		}
		return GroupIterator(p);
	}

	GroupIterator make_end_iterator(std::string const& familyName) {
		boost::shared_ptr<AttributeRangeIterator<T> > p;
		Simulation::const_iterator simIter = sim.find(familyName);
		if(simIter == sim.end())
			p.reset(new AttributeRangeIterator<T>(0, 0, 0, 0, 0));
		else {
			AttributeMap::const_iterator attrIter = simIter->second.attributes.find(attributeName);
			TypedArray const& array = attrIter->second;
			p.reset(new AttributeRangeIterator<T>(array.length, minValue, maxValue, array.getArray(Type2Type<T>()), array.length));
		}
		return GroupIterator(p);
	}
};

boost::shared_ptr<Group> make_AttributeRangeGroup(Simulation const& sim, std::string const& attributeName, double minValue, double maxValue) {
	boost::shared_ptr<Group> p;
	for(Simulation::const_iterator simIter = sim.begin(); simIter != sim.end(); ++simIter) {
		AttributeMap::const_iterator attrIter = simIter->second.attributes.find(attributeName);
		if(attrIter != simIter->second.attributes.end()) {
			//found a family with the attribute
			TypedArray const& arr = attrIter->second;
			if(arr.dimensions == 1) {
				switch(arr.code) {
					case int8:
						p.reset(new AttributeRangeGroup<Code2Type<int8>::type>(sim, attributeName, static_cast<Code2Type<int8>::type>(minValue), static_cast<Code2Type<int8>::type>(maxValue)));
						break;
					case uint8:
						p.reset(new AttributeRangeGroup<Code2Type<uint8>::type>(sim, attributeName, static_cast<Code2Type<uint8>::type>(minValue), static_cast<Code2Type<uint8>::type>(maxValue)));
						break;
					case int16:
						p.reset(new AttributeRangeGroup<Code2Type<int16>::type>(sim, attributeName, static_cast<Code2Type<int16>::type>(minValue), static_cast<Code2Type<int16>::type>(maxValue)));
						break;
					case uint16:
						p.reset(new AttributeRangeGroup<Code2Type<uint16>::type>(sim, attributeName, static_cast<Code2Type<uint16>::type>(minValue), static_cast<Code2Type<uint16>::type>(maxValue)));
						break;
					case int32:
						p.reset(new AttributeRangeGroup<Code2Type<int32>::type>(sim, attributeName, static_cast<Code2Type<int32>::type>(minValue), static_cast<Code2Type<int32>::type>(maxValue)));
						break;
					case uint32:
						p.reset(new AttributeRangeGroup<Code2Type<uint32>::type>(sim, attributeName, static_cast<Code2Type<uint32>::type>(minValue), static_cast<Code2Type<uint32>::type>(maxValue)));
						break;
					case int64:
						p.reset(new AttributeRangeGroup<Code2Type<int64>::type>(sim, attributeName, static_cast<Code2Type<int64>::type>(minValue), static_cast<Code2Type<int64>::type>(maxValue)));
						break;
					case uint64:
						p.reset(new AttributeRangeGroup<Code2Type<uint64>::type>(sim, attributeName, static_cast<Code2Type<uint64>::type>(minValue), static_cast<Code2Type<uint64>::type>(maxValue)));
						break;
					case float32:
						p.reset(new AttributeRangeGroup<Code2Type<float32>::type>(sim, attributeName, static_cast<Code2Type<float32>::type>(minValue), static_cast<Code2Type<float32>::type>(maxValue)));
						break;
					case float64:
						p.reset(new AttributeRangeGroup<Code2Type<float64>::type>(sim, attributeName, static_cast<Code2Type<float64>::type>(minValue), static_cast<Code2Type<float64>::type>(maxValue)));
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
		
	FamilyGroup(Simulation const& s, std::string const& familyName) : sim(s) {
		families.push_back(familyName);
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

//Static groups

//Tree-based groups

} //close namespace SimulationHandling

#endif //GROUP_H__7h237h327yh4th7y378w37aw7ho4wvta8