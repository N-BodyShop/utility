/** \file SPH_Kernel.h
 \author Graeme Lufkin (gwl@u.washington.edu)
 \date Created Summer 2002
 \version 1.0
 */

#ifndef SPH_KERNEL_H
#define SPH_KERNEL_H

namespace SPH {

class Kernel {
public:
	/** Evaluate the kernel for a given distance r with a given smoothing length h.
	 The distance r is in the range [0, 2h]
	 */
	virtual double evaluate(double r, double h) const = 0;
	
	/** Evaluate the gradient (almost) of the kernel for a given distance r with a given smoothing length h.
	 The distance r is in the range [0, 2h]
	 To get the full gradient, multiply the result of this function by the vector whose magnitude is r.
	 */
	virtual double evaluateGradient(double r, double h) const = 0;
};

/** The standard cubic kernel used in SPH calculations.
 */
class SplineKernel : public Kernel {
	static const double spline_PI = 3.14159265358979323846;
public:

	inline double evaluate(double r, double h) const {
		double q = r / h;
		if(q < 1)
			return (1 - 1.5 * q * q + 0.75 * q * q * q) / spline_PI / h / h / h;
		else
			return 0.25 * (2 - q) * (2 - q) * (2 - q) / spline_PI / h / h / h;
	}
	
	inline double evaluateGradient(double r, double h) const {
		double q = r / h;
		if(q < 1)
			return (0.75 * q - 1) * 3 / spline_PI / h / h / h / h / h;
		else
			return (-0.25 * q - 1 / q + 1) * 3 / spline_PI / h / h / h / h / h;
	}
};

} //close namespace SPH

#endif //SPH_KERNEL_H
