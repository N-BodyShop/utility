/** \file OrientedBox.h
 \author Graeme Lufkin (gwl@u.washington.edu)
 \date Created October 22, 2001
 \version 1.0
 \todo Fully document this class
 */

#ifndef ORIENTEDBOX_H
#define ORIENTEDBOX_H

#include <iostream>
#include <sstream>
#include <string>

#include "Vector3D.h"
#include "Shape.h"
#include "Sphere.h"
#include "Box.h"
#include "PeriodicBoundaryConditions.h"

/// A box in three dimensions whose axes are aligned with the coordinate axes
template <class T>
class OrientedBox : public Shape<T> {
private:

public:
	/// The corner with the minimum \c x, \c y, and \c z values
	Vector3D<T> lesser_corner;
	/// The corner with the maximum \c x, \c y, and \c z values
	Vector3D<T> greater_corner;
	
	OrientedBox(const Vector3D<T>& corner1 = Vector3D<T>(-0.5, -0.5, -0.5), const Vector3D<T>& corner2 = Vector3D<T>(0.5, 0.5, 0.5)) {
		if(corner1.x < corner2.x) {
			lesser_corner = Vector3D<T>(corner1);
			if(corner1.y > corner2.y || corner1.z > corner2.z) //malformed box!
				greater_corner = lesser_corner;
			else
				greater_corner = Vector3D<T>(corner2);
		} else {
			lesser_corner = Vector3D<T>(corner2);
			if(corner1.y < corner2.y || corner1.z < corner2.z) //malformed box!
				greater_corner = lesser_corner;
			else
				greater_corner = Vector3D<T>(corner1);
		}
		
	}
	
	OrientedBox(const Box<T>& b) {
		lesser_corner = greater_corner = b.vertices[0];
		for(int i = 0; i < 8; i++) {
			if(b.vertices[i].x < lesser_corner.x)
				lesser_corner.x = b.vertices[i].x;
			if(b.vertices[i].x > greater_corner.x)
				greater_corner.x = b.vertices[i].x;
			if(b.vertices[i].y < lesser_corner.y)
				lesser_corner.y = b.vertices[i].y;
			if(b.vertices[i].y > greater_corner.y)
				greater_corner.y = b.vertices[i].y;
			if(b.vertices[i].z < lesser_corner.z)
				lesser_corner.z = b.vertices[i].z;
			if(b.vertices[i].z > greater_corner.z)
				greater_corner.z = b.vertices[i].z;
		}
	}
	
	~OrientedBox() { }

	inline bool contains(const Vector3D<T>& point) const {
		return point.x >= lesser_corner.x && point.x <= greater_corner.x
				&& point.y >= lesser_corner.y && point.y <= greater_corner.y
				&& point.z >= lesser_corner.z && point.z <= greater_corner.z;
	}
	
	inline Vector3D<T> center() const {
		return (lesser_corner + greater_corner) / 2.0;
	}
	
	inline void grow(const Vector3D<T>& point) {
		if(point.x < lesser_corner.x)
			lesser_corner.x = point.x;
		if(point.y < lesser_corner.y)
			lesser_corner.y = point.y;
		if(point.z < lesser_corner.z)
			lesser_corner.z = point.z;
		
		if(point.x > greater_corner.x)
			greater_corner.x = point.x;
		if(point.y > greater_corner.y)
			greater_corner.y = point.y;
		if(point.z > greater_corner.z)
			greater_corner.z = point.z;		
	}
	
	inline T volume() const {
		return (greater_corner.x - lesser_corner.x) * (greater_corner.y - lesser_corner.y) * (greater_corner.z - lesser_corner.z);
	}
	
	template <class T2>
	inline bool encloses(const OrientedBox<T2>& b) const {
		return lesser_corner.x <= b.lesser_corner.x
				&& lesser_corner.y <= b.lesser_corner.y
				&& lesser_corner.z <= b.lesser_corner.z
				&& greater_corner.x >= b.greater_corner.x
				&& greater_corner.y >= b.greater_corner.y
				&& greater_corner.z >= b.greater_corner.z;
	}
	
	template <class T2>
	inline bool intersects(const OrientedBox<T2>& b) const {
		return !(lesser_corner.x > b.greater_corner.x || greater_corner.x < b.lesser_corner.x
				|| lesser_corner.y > b.greater_corner.y || greater_corner.y < b.lesser_corner.y
				|| lesser_corner.z > b.greater_corner.z || greater_corner.z < b.lesser_corner.z);
	}

	/// Determine if a sphere and a box intersect
	template <class T2>
	inline bool intersects(const Sphere<T2>& s) const {
		T dsq = 0;
		T rsq = s.radius * s.radius;
		T delta;
		if((delta = lesser_corner.x - s.origin.x) > 0)
			dsq += delta * delta;
		else if((delta = s.origin.x - greater_corner.x) > 0)
			dsq += delta * delta;
		if(rsq < dsq)
			return false;
		if((delta = lesser_corner.y - s.origin.y) > 0)
			dsq += delta * delta;
		else if((delta = s.origin.y - greater_corner.y) > 0)
			dsq += delta * delta;
		if(rsq < dsq)
			return false;
		if((delta = lesser_corner.z - s.origin.z) > 0)
			dsq += delta * delta;
		else if((delta = s.origin.z - greater_corner.z) > 0)
			dsq += delta * delta;
		return (dsq <= s.radius * s.radius);
	}

	/// Determine if a sphere and a box intersect, using periodic boundary conditions
	template <class T2>
	inline bool intersects(const PeriodicBoundaryConditions<T2>& pbc, const Sphere<T2>& s) const {
		T dsq = 0;
		T rsq = s.radius * s.radius;
		T delta, delta_ghost;
		if((delta = lesser_corner.x - s.origin.x) > 0) {
			if(pbc.xPeriod > 0) {
				delta_ghost = s.origin.x + pbc.xPeriod - greater_corner.x;
				delta = (delta < delta_ghost ? delta : delta_ghost);
			}
			dsq += delta * delta;
		} else if((delta = s.origin.x - greater_corner.x) > 0) {
			if(pbc.xPeriod > 0) {
				delta_ghost = lesser_corner.x - s.origin.x + pbc.xPeriod;
				delta = (delta < delta_ghost ? delta : delta_ghost);
			}
			dsq += delta * delta;
		}
		if(rsq < dsq)
			return false;
		if((delta = lesser_corner.y - s.origin.y) > 0) {
			if(pbc.yPeriod > 0) {
				delta_ghost = s.origin.y + pbc.yPeriod - greater_corner.y;
				delta = (delta < delta_ghost ? delta : delta_ghost);
			}
			dsq += delta * delta;
		} else if((delta = s.origin.y - greater_corner.y) > 0) {
			if(pbc.yPeriod > 0) {
				delta_ghost = lesser_corner.y - s.origin.y + pbc.yPeriod;
				delta = (delta < delta_ghost ? delta : delta_ghost);
			}
			dsq += delta * delta;
		}
		if(rsq < dsq)
			return false;
		if((delta = lesser_corner.z - s.origin.z) > 0) {
			if(pbc.zPeriod > 0) {
				delta_ghost = s.origin.z + pbc.zPeriod - greater_corner.z;
				delta = (delta < delta_ghost ? delta : delta_ghost);
			}
			dsq += delta * delta;
		} else if((delta = s.origin.z - greater_corner.z) > 0) {
			if(pbc.zPeriod > 0) {
				delta_ghost = lesser_corner.z - s.origin.z + pbc.zPeriod;
				delta = (delta < delta_ghost ? delta : delta_ghost);
			}
			dsq += delta * delta;
		}
		return (dsq <= s.radius * s.radius);
	}
	
	/// Output operator, used for formatted display
	friend std::ostream& operator<< (std::ostream& os, const OrientedBox<T>& b) {
		os << '{' << b.lesser_corner << ',' << b.greater_corner << '}';
		return os;
	}

	/** Make string containing the tipsy command to create this box.
	 \param i The number of the box to make
	 \return a string that can be given to tipsy
	 */
	std::string tipsyCommand(int i) const {
		std::ostringstream oss;
		oss << "setbox " << i << " ";
		oss << (greater_corner.x / 2 + lesser_corner.x / 2) << " "
				<< (greater_corner.y / 2 + lesser_corner.y / 2) << " "
				<< (greater_corner.z / 2 + lesser_corner.z / 2) << " "
				<< (greater_corner.x / 2 - lesser_corner.x / 2) << " "
				<< (greater_corner.y / 2 - lesser_corner.y / 2) << " "
				<< (greater_corner.z / 2 - lesser_corner.z / 2);
		return oss.str();
	}
};

#endif //ORIENTEDBOX_H
