//export_Vector3D.cpp

#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/operators.hpp>

#include "Vector3D.h"

using namespace boost::python;

typedef Vector3D<float> Vec;

std::string bare(const Vec& v) {
	std::ostringstream oss;
	oss << v.x << " " << v.y << " " << v.z;
	return oss.str();
}

void export_Vector3D() {
	
	class_<Vec>("Vector3D", "A Cartesian vector in three dimensions", init<float, float, float>())
		.def(init<>())
		.def(init<const Vec&>())
		.def_readwrite("x", &Vec::x)
		.def_readwrite("y", &Vec::y)
		.def_readwrite("z", &Vec::z)
		.def("length", &Vec::length)
		.def("lengthSquared", &Vec::lengthSquared)
		.def("normalize", &Vec::normalize, return_value_policy<reference_existing_object>())
		.def(self == self)
		.def(self != self)
		.def(self + self)
		.def(self - self)
		.def(self += self)
		.def(self -= self)
		.def(self * float())
		.def(self *= float())
		.def(self / float())
		.def(self /= float())
		.def(-self)
		.def(self * self)
		.def(self *= self)
		.def(self / self)
		.def(self /= self)
		.def(float() * self)
		.def(str(self))
		.def("bare", bare, "A string containing the components of this vector")
		;
		
	def("dot", static_cast<float (*)(const Vec&, const Vec&)>(dot), "The dot product of two vectors");
	def("costheta",  static_cast<float (*)(const Vec&, const Vec&)>(costheta), "The cosine of the angle between two vectors");
	def("cross",  static_cast<Vec (*)(const Vec&, const Vec&)>(cross), "The cross product of two vectors");
	
	def("fromSpherical",  static_cast<Vec (*)(const float, const float, const float)>(fromSpherical), (arg("r"), arg("theta"), arg("phi")), "Construct a Cartesian vector from spherical components");
	def("fromCylindrical",  static_cast<Vec (*)(const float, const float, const float)>(fromCylindrical), (arg("r"), arg("theta"), arg("z")), "Construct a Cartesian vector from cylindrical components");
	
	class_<RotationMatrix<float> >("RotationMatrix", init<const float, const float, const float>())
		.def("rotate", &RotationMatrix<float>::rotate)
		;

}
