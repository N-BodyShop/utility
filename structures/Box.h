/** \file Box.h
 \author Graeme Lufkin (gwl@u.washington.edu)
 \date Created December 3, 2001
 \version 1.0
 \todo Implement this class!
 */

#ifndef BOX_H
#define BOX_H

#include "Vector3D.h"
#include "Shape.h"

template <class T>
class Face {
private:
	Vector3D<T> center;
	T width;
	T height;
	Vector3D<T> norm;
public:
	Face() { }

	~Face() { }
	
	Vector3D<T> normal() const {
		return norm;
	}
};

template <class T>
class Box : public Shape<T> {
private:
	Face<T> faces[6];
	T length, width, height;
	
public:
	Vector3D<T> vertices[8];

	Box() { }

	Box(const Vector3D<T>& corner0, const Vector3D<T>& corner1,
		const Vector3D<T>& corner2, const Vector3D<T>& corner3) {
	}
	
	~Box() { }
	
	bool contains(const Vector3D<T>& point) const {
		return 1;
	}
	
	void grow(const Vector3D<T>& point) {
		
	}
	
	T volume() const {
		return length * width * height;
	}
};

#endif //BOX_H
