/** \file Shape.h
 \author Graeme Lufkin (gwl@u.washington.edu)
 \date Created November 8, 2001
 \version 1.0
 */

#ifndef SHAPE_H
#define SHAPE_H

#include "Vector3D.h"

/// The base class defining a geometric shape in three dimensions
template <class T>
class Shape {
public:
	Shape() { }
	virtual ~Shape() { };
	
	/** Does this shape contain (enclose) a particular point?
	 \param point The point of interest
	 \return \c true if the shape contains the point, \c false otherwise
	 */
	virtual bool contains(const Vector3D<T>& point) const = 0;
	
	/// Make the shape grow larger (if necessary) to accomodate this point
	virtual void grow(const Vector3D<T>& point) = 0;
	
	/// What is the volume enclosed by this shape?
	virtual T volume() const = 0;
};

#endif //SHAPE_H
