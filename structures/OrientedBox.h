/** @file OrientedBox.h
 @author Graeme Lufkin (gwl@u.washington.edu)
 @date Created October 22, 2001
 @version 1.0
 @todo Fully document this class
 */

#ifndef ORIENTEDBOX_H
#define ORIENTEDBOX_H

#include <iostream>
#include <string>
#include <sstream>

#include "Shape.h"
#include "Vector3D.h"
#include "Sphere.h"

/// A box in three dimensions whose axes are aligned with the coordinate axes
template <typename T = double>
class OrientedBox : public Shape<T> {
private:

public:
	/// The corner with the minimum \c x, \c y, and \c z values
	Vector3D<T> lesser_corner;
	/// The corner with the maximum \c x, \c y, and \c z values
	Vector3D<T> greater_corner;
	
	OrientedBox() {
		lesser_corner = Vector3D<T>(HUGE_VAL, HUGE_VAL, HUGE_VAL);
		greater_corner = Vector3D<T>(-HUGE_VAL, -HUGE_VAL, -HUGE_VAL);
	}
	
	explicit OrientedBox(const Vector3D<T>& corner1, const Vector3D<T>& corner2) {
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
	
	OrientedBox(const OrientedBox<T>& b) : lesser_corner(b.lesser_corner), greater_corner(b.greater_corner) { }
	
	~OrientedBox() { }
	
	OrientedBox<T>& operator=(const OrientedBox<T>& b) {
		lesser_corner = b.lesser_corner;
		greater_corner = b.greater_corner;
		return *this;
	}
	
	template <class T2>
	operator OrientedBox<T2> () const {
		return OrientedBox<T2>(static_cast<Vector3D<T2> >(lesser_corner), static_cast<Vector3D<T2> >(greater_corner));
	}
	
	inline bool contains(const Vector3D<T>& point) const {
		return point.x >= lesser_corner.x && point.x <= greater_corner.x
				&& point.y >= lesser_corner.y && point.y <= greater_corner.y
				&& point.z >= lesser_corner.z && point.z <= greater_corner.z;
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
	
	inline Vector3D<T> center() const {
		return (lesser_corner + greater_corner) / 2.0;
	}
	
	inline Vector3D<T> size() const {
		return greater_corner - lesser_corner;
	}
	
	inline OrientedBox<T>& shift(const Vector3D<T>& v) {
		lesser_corner += v;
		greater_corner += v;
		return *this;
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

#ifdef CHARM
#include "pup.h"

template <typename T>
inline void operator|(PUP::er& p, OrientedBox<T>& b) {
	p | b.lesser_corner;
	p | b.greater_corner;
}

#endif //CHARM

#endif //ORIENTEDBOX_H
