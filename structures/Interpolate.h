/** \file Interpolate.h
 This file defines classes that represent interpolated functions. 
 \author Graeme Lufkin (gwl@u.washington.edu)
 \date Created May 1, 2000
 \version 2.0
 \todo Fully document this file.
 */

#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include <vector>
#include <algorithm>

using std::vector;
using std::copy;

/**
 * This class does linear interpolation on a data set with function values 
 * specified initially at a number of values of the independent variable.
 * The independent variable values need not be equally spaced.
 * After initialization, calls to this object with a value of the independent
 * variable will return the linear interpolation of the value of the function
 * at the requested point, based on the surrounding values.  If you ask for
 * a value outside the initially specified range, the end value will be returned.
 */
template <class T>
class LinearInterpolator {
private:
	
	bool ready;
	vector<T> indeps;
	vector<T> deps;
	vector<T> coefs;
	unsigned int size;
	mutable int klo, khi;
	
public:
	
	/// Default constructor, is non-working
	LinearInterpolator() {
		ready = false;
	}
	
	//create a linear interpolator object, given iterator pairs for the independent
	//variables and the dependent variables (the function values at the specified points).
	template <typename InputIterator, typename InputIterator2>
	LinearInterpolator(InputIterator beginIndep, InputIterator endIndep, 
			InputIterator2 beginDep, InputIterator2 endDep) : indeps(beginIndep, endIndep), 
			deps(beginDep, endDep) {

		size = indeps.size();
		if(size != deps.size() || size == 0) {
			ready = false;
			return;
		}
		
		coefs.assign(size - 1, 0);
		
		//calculate the coefficients to be used when interpolating
		for(unsigned int i = 0; i < size - 1; i++)
			coefs[i] = (deps[i + 1] - deps[i]) / (indeps[i + 1] - indeps[i]);
		
		klo = 0;
		khi = size - 1;		
		
		ready = true;
	}
	
	/// The linear interpolator object can be called just like a regular function
	T operator()(T x) const {
		if(x < indeps[klo]) { //lower bracket than last time
			if(x < indeps[0])
				return deps[0];
			else if(indeps[klo - 1] < x) { //it's the next one down
				klo--;
				khi--;
			} else //set low location as low as possible
				klo = 0;
		} else if(indeps[khi] < x) { //greater bracket than last time
			if(indeps[size - 1] < x)
				return deps[size - 1];
			else if(x < indeps[khi + 1]) { //it's the next one up
				khi++;
				klo++;
			} else //set high location as high as possible
				khi = size - 1;
		} //else, it's the same bracket as last time, don't change bounds
		
		int k;
		while(khi - klo > 1) { //do binary search to find bracket
			k = (khi + klo) / 2;
			if(x < indeps[k])
				khi = k;
			else
				klo = k;
		}
		
		//use the coefficients to give f(x) = a * x + b
		return coefs[klo] * (x - indeps[klo]) + deps[klo];
	}
	
	bool isReady() const {
		return ready;
	}
};

template <class T>
class SplineDerivative;

/**
 * This class does cubic spline interpolation on a data set with function values 
 * specified initially at a number of values of the independent variable.
 * It functions identically to the linear interpolator.  The first derivative
 * at the first and last points can be specified.  If they are not, then the
 * "natural" spline is used, where the second derivative is zero at the end
 * points.  If you ask for a value outside the initially specified range,
 * the linear extrapolation is returned, using the value of the first 
 * derivative at the end point.
 */
template <class T>
class SplineInterpolator {
	friend class SplineDerivative<T>;
private:
	
	bool ready;
	vector<T> indeps;
	vector<T> deps;
	vector<T> secondDerivs;
	T beginFirstDeriv;
	T endFirstDeriv;
	unsigned int size;
	mutable int klo, khi;
	
public:
	
	//default constructor, is non-working
	SplineInterpolator() {
		ready = false;
	}
		
	/* Create a spline interpolator object.
	 * This version takes specific values for the first derivatives at the end points.
	 */
	template <typename InputIterator, typename InputIterator2>
	SplineInterpolator(InputIterator beginIndep, InputIterator endIndep,
			InputIterator2 beginDep, InputIterator2 endDep, T beginPrime,
			T endPrime) : indeps(beginIndep, endIndep), deps(beginDep, endDep),
			beginFirstDeriv(beginPrime), endFirstDeriv(endPrime) {
		
		size = indeps.size();
		if(size != deps.size() || size == 0) {
			ready = false;
			return;
		}
		
		secondDerivs.assign(size, 0);
		
		vector<T> u(size - 1);
		T sig, p, qn, un;
		
		secondDerivs[0] = -0.5;
		u[0] = (3.0 / (indeps[1] - indeps[0])) * ((deps[1] - deps[0]) / (indeps[1] - indeps[0]) - beginFirstDeriv);

		for(int i = 1; i < size - 1; i++) {
			sig = (indeps[i] - indeps[i - 1]) / (indeps[i + 1] - indeps[i - 1]);
			p = sig * secondDerivs[i - 1] + 2.0;
			secondDerivs[i] = (sig - 1.0) / p;
			u[i] = (deps[i + 1] - deps[i]) / (indeps[i + 1] - indeps[i]) - (deps[i] - deps[i - 1]) / (indeps[i] - indeps[i - 1]);
			u[i] = (6.0 * u[i] / (indeps[i + 1] - indeps[i - 1]) - sig * u[i - 1]) / p;
		}
		
		qn = 0.5;
		un = (3.0 / (indeps[size - 1] - indeps[size - 2])) * (endFirstDeriv - (deps[size - 1] - deps[size - 2]) / (indeps[size - 1] - indeps[size - 2]));

		secondDerivs[size - 1] = (un - qn * u[size - 2]) / (qn * secondDerivs[size - 2] + 1.0);
		for(int i = size - 2; i >= 0; i--)
			secondDerivs[i] = secondDerivs[i] * secondDerivs[i + 1] + u[i];
				
		klo = 0;
		khi = size - 1;		
		
		ready = true;
	}
	
	/* Create a spline interpolator object.
	 * This version creates a "natural" spline, with the second derivative
	 * equal to zero at the end points.
	 */
	template <typename InputIterator, typename InputIterator2>
	SplineInterpolator(InputIterator beginIndep, InputIterator endIndep,
			InputIterator2 beginDep, InputIterator2 endDep) : indeps(beginIndep, endIndep),
			deps(beginDep, endDep) {
		
		size = indeps.size();
		if(size != deps.size() || size == 0) {
			ready = false;
			return;
		}

		//calculate the value of the first derivative at the end points
		beginFirstDeriv = (deps[1] - deps[0]) / (indeps[1] - indeps[0]);
		endFirstDeriv = (deps[size - 1] - deps[size - 2]) / (indeps[size - 1] - indeps[size - 2]);
		
		secondDerivs.assign(size, 0);
		
		vector<T> u(size - 1);
		T sig, p;
		
		secondDerivs[0] = 0;
		u[0] = 0;
		
		for(unsigned int i = 1; i < size - 1; i++) {
			sig = (indeps[i] - indeps[i - 1]) / (indeps[i + 1] - indeps[i - 1]);
			p = sig * secondDerivs[i - 1] + 2.0;
			secondDerivs[i] = (sig - 1.0) / p;
			u[i] = (deps[i + 1] - deps[i]) / (indeps[i + 1] - indeps[i]) - (deps[i] - deps[i - 1]) / (indeps[i] - indeps[i - 1]);
			u[i] = (6.0 * u[i] / (indeps[i + 1] - indeps[i - 1]) - sig * u[i - 1]) / p;
		}

		//copy(u.begin(), u.end(), ostream_iterator<double>(cerr, " "));
		//cerr << endl;
		
		secondDerivs[size - 1] = 0;
		for(int i = size - 2; i >= 0; i--)
			secondDerivs[i] = secondDerivs[i] * secondDerivs[i + 1] + u[i];
		
		klo = 0;
		khi = size - 1;		
		
		ready = true;
	}
	
	//the spline interpolator object can be called just like a regular function
	T operator()(T x) const {
		if(x < indeps[klo]) { //lower bracket than last time
			if(x < indeps[0]) //before the beginning
				return deps[0] - beginFirstDeriv * (indeps[0] - x);
			else if(indeps[klo - 1] < x) { //it's the next one down
				klo--;
				khi--;
			} else //set low location as low as possible
				klo = 0;
		} else if(indeps[khi] < x) { //greater bracket than last time
			if(indeps[size - 1] < x) //beyond the end
				return deps[size - 1] + endFirstDeriv * (x - indeps[size - 1]);
			else if(x < indeps[khi + 1]) { //it's the next one up
				khi++;
				klo++;
			} else //set high location as high as possible
				khi = size - 1;
		} //else, it's the same bracket as last time, don't change bounds

		int k;
		while(khi - klo > 1) { //do binary search to find bracket
			k = (khi + klo) / 2;
			if(x < indeps[k])
				khi = k;
			else
				klo = k;
		}
		
		//return the interpolation of the function
		T h = indeps[khi] - indeps[klo];
		T A = (indeps[khi] - x) / h;
		T B = (x - indeps[klo]) / h;
		return A * deps[klo] + B * deps[khi] 
				+ h * h * (A * A * A - A) * secondDerivs[klo] / 6.0
				+ h * h * (B * B * B - B) * secondDerivs[khi] / 6.0;
	}
	
	bool isReady() const {
		return ready;
	}
};

/**
 * This class represents the derivative of a cubic-spline interpolated function.
 */
template <class T>
class SplineDerivative {
private:
	
	bool ready;
	vector<T> indeps;
	vector<T> deps;
	vector<T> secondDerivs;
	T beginFirstDeriv;
	T endFirstDeriv;
	unsigned int size;
	mutable int klo, khi;
	
public:
	
	//default constructor, is non-working
	SplineDerivative() {
		ready = false;
	}
	
	SplineDerivative(const SplineInterpolator<T>& si) : indeps(si.indeps), deps(si.deps), 
			secondDerivs(si.secondDerivs), beginFirstDeriv(si.beginFirstDeriv),
			endFirstDeriv(si.endFirstDeriv) {
		
		size = indeps.size();
		if(size != deps.size()) {
			ready = false;
			return;
		}
		
		klo = 0;
		khi = size - 1;		
		
		ready = true;
	}
	
	//the spline derivative object can be called just like a regular function
	T operator()(T x) const {
		if(x < indeps[klo]) { //lower bracket than last time
			if(x < indeps[0]) //before the beginning
				return beginFirstDeriv;
			else if(indeps[klo - 1] < x) { //it's the next one down
				klo--;
				khi--;
			} else //set low location as low as possible
				klo = 0;
		} else if(indeps[khi] < x) { //greater bracket than last time
			if(indeps[size - 1] < x) //beyond the end
				return endFirstDeriv;
			else if(x < indeps[khi + 1]) { //it's the next one up
				khi++;
				klo++;
			} else //set high location as high as possible
				khi = size - 1;
		} //else, it's the same bracket as last time, don't change bounds

		int k;
		while(khi - klo > 1) { //do binary search to find bracket
			k = (khi + klo) / 2;
			if(x < indeps[k])
				khi = k;
			else
				klo = k;
		}
		
		//return the derivative of the interpolation of the function
		
		T h = indeps[khi] - indeps[klo];
		T A = (indeps[khi] - x) / h;
		T B = (x - indeps[klo]) / h;
		return (deps[khi] - deps[klo]) / h 
				- (3 * A * A * A - 1) * h * secondDerivs[klo] / 6
				+ (3 * B * B * B - 1) * h * secondDerivs[khi] / 6;
	}

	bool isReady() const {
		return ready;
	}
};

#endif //INTERPOLATE_H
