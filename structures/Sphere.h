/** \file Shape.h
 \author Graeme Lufkin (gwl@u.washington.edu)
 \date Created November 2, 2001
 \version 1.0
 */

#ifndef SPHERE_H
#define SPHERE_H

#include <iostream>

#include "Vector3D.h"
#include "PeriodicBoundaryConditions.h"

template <typename T>
class OrientedBox;

/// A class representing a sphere in three dimensions
template <typename T = double>
class Sphere {
public:
	/// The origin of this sphere
	Vector3D<T> origin;
	/// The radius of this sphere
	T radius;

	Sphere(const Vector3D<T>& o = Vector3D<T>(), const T r = 1) : origin(o), radius(r) { }
	
	Sphere(const Sphere<T>& s) : origin(s.origin), radius(s.radius) { }
	
	~Sphere() { }
	
	Sphere<T>& operator=(const Sphere<T>& s) {
		origin = s.origin;
		radius = s.radius;
		return *this;
	}
	
	/// A sphere contains a point if the distance between the origin and the point is less than the radius
	inline bool contains(const Vector3D<T>& point) const {
		return (origin - point).length() < radius;
	}

	/// Does this sphere contain a point, given periodic boundary conditions
	inline bool contains(const PeriodicBoundaryConditions<T>& pbc, const Vector3D<T>& point) const {
		T dsq = 0;
		T rsq = radius * radius;
		T delta, delta_ghost;
		delta = origin.x - point.x;
		if(pbc.xPeriod > 0) {
			if(delta > 0)
				delta_ghost = point.x + pbc.xPeriod - origin.x;
			else {
				delta *= -1;
				delta_ghost = origin.x - point.x + pbc.xPeriod;
			}
			delta = (delta < delta_ghost ? delta : delta_ghost);
		}
		dsq += delta * delta;
		if(rsq < dsq)
			return false;
		delta = origin.y - point.y;
		if(pbc.yPeriod > 0) {
			if(delta > 0)
				delta_ghost = point.y + pbc.yPeriod - origin.y;
			else {
				delta *= -1;
				delta_ghost = origin.y - point.y + pbc.yPeriod;
			}
			delta = (delta < delta_ghost ? delta : delta_ghost);
		}
		dsq += delta * delta;
		if(rsq < dsq)
			return false;
		delta = origin.z - point.z;
		if(pbc.zPeriod > 0) {
			if(delta > 0)
				delta_ghost = point.z + pbc.zPeriod - origin.z;
			else {
				delta *= -1;
				delta_ghost = origin.z - point.z + pbc.zPeriod;
			}
			delta = (delta < delta_ghost ? delta : delta_ghost);
		}
		dsq += delta * delta;
		return (dsq < rsq);
	}

	/// Growing a sphere keeps the origin fixed and increases the radius
	inline void grow(const Vector3D<T>& point) {
		T points_radius = (origin - point).length();
		if(points_radius > radius)
			radius = points_radius;
	}
	
	/// The volume of a sphere is \f$ \frac{4 \pi}{3} r^3 \f$
	inline T volume() const {
		return 4.0 * 3.14159265358979323846 / 3.0 * radius * radius * radius;
	}
	
	/// Two spheres intersect if the distance betweeen their origins is larger than the sum of their radii
	inline bool intersects(const Sphere<T>& s) const {
		return (origin - s.origin).length() <= (radius + s.radius);
	}
	
	/// Output operator, used for formatted display
	friend std::ostream& operator<<(std::ostream& os, const Sphere<T>& s) {
		os << '{' << s.origin << ", " << s.radius << '}';
		return os;
	}

};

typedef Sphere<double> dSphere;

#endif //SPHERE_H
