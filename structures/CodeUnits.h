//CodeUnits.h

#ifndef CODEUNITS_H__hsdf07j253v89724tvuhi23r2t3vuj23v
#define CODEUNITS_H__hsdf07j253v89724tvuhi23r2t3vuj23v

#include <cmath>

/** This class encapsulates the conversion factors from a set of simulation units to CGS. 
 To use the conversions, instantiate a \c CodeUnits object with the unit of length in
 centimeters and the unit of mass in grams.  By requiring that the gravitational
 constant G be equal to 1, this specifies the unit of time.  Several useful additional
 conversion factors are derived from the three fundamental conversions.  Temperatures
 are always measured in Kelvin.
 To convert a quantity from code units to CGS, multiply it by the conversion factor
 for that unit.  To go from CGS to code units, divide the quantity by the factor.
*/
class CodeUnits {
public:
	
	//Some physical constants in CGS units
	static const double BoltzmannConstantCGS = 1.38065E-16; //k_b in grams cm^2 per Kelvin per second^2
	static const double StefanBoltzmannConstantCGS = 5.6704E-5;  // \sigma in grams per Kelvin^4 per second^3
	static const double GravitationalConstantCGS = 6.673E-8;
	static const double HydrogenMassCGS = 1.67339E-24; //m_H
	static const double GasConstantCGS = 8.2506188037E7; //k_B / m_H
	static const double AU = 1.4959787066E13; //in cm
	static const double KiloParsec = 3.0857E21; //in cm
	static const double EarthMass = 5.9742E27; //in grams
	static const double SolarMass = 1.9891E33; //in grams
	static const double Year = 3.1536E7; //in seconds
	
	//force G == 1 in code units
	static const double GravitationalConstant = 1.0;

	//Conversions for the unit of distance to commonly used units
	double distanceUnitInCentimeters;
	double distanceUnitInAU;
	double distanceUnitInKpc;
	
	//Conversions for the unit of mass to commonly used units
	double massUnitInGrams;
	double massUnitInEarthMasses;
	double massUnitInSolarMasses;
	
	//Conversions for the unit of time to commonly used units
	double timeUnitInSeconds;
	double timeUnitInYears;
	
	//Conversions for various quantities, built out of the fundamental conversions
	double densityUnitCGS;
	double forceUnitsCGS;
	double pressureUnitCGS;
	double opacityUnitCGS;
	
	// Values of physical constants in code units
	double BoltzmannConstant;
	double StefanBoltzmannConstant;
	double HydrogenMass;
	double GasConstant;
	
	/// Define the system of units by specifying the unit of distance in centimeters and the unit of mass in grams (by requiring that G = 1, this specifies the unit of time).
	CodeUnits(const double distanceUnit = 0, const double massUnit = 0) : 
			distanceUnitInCentimeters(distanceUnit), 
			massUnitInGrams(massUnit), 
			timeUnitInSeconds(sqrt(distanceUnit * distanceUnit * distanceUnit / massUnit / GravitationalConstantCGS))
			{
		
		distanceUnitInAU = distanceUnitInCentimeters / AU;
		distanceUnitInKpc = distanceUnitInCentimeters / KiloParsec;
		massUnitInEarthMasses = massUnitInGrams / EarthMass;
		massUnitInSolarMasses = massUnitInGrams / SolarMass;
		timeUnitInYears = timeUnitInSeconds / Year;
		
		densityUnitCGS = massUnitInGrams / distanceUnitInCentimeters / distanceUnitInCentimeters / distanceUnitInCentimeters;
		forceUnitsCGS = massUnitInGrams * distanceUnitInCentimeters / timeUnitInSeconds / timeUnitInSeconds;
		pressureUnitCGS = massUnitInGrams / distanceUnitInCentimeters / timeUnitInSeconds / timeUnitInSeconds;
		opacityUnitCGS = distanceUnitInCentimeters * distanceUnitInCentimeters / massUnitInGrams;
		
		BoltzmannConstant = BoltzmannConstantCGS / massUnitInGrams / distanceUnitInCentimeters / distanceUnitInCentimeters * timeUnitInSeconds * timeUnitInSeconds;
		StefanBoltzmannConstant = StefanBoltzmannConstantCGS / massUnitInGrams * timeUnitInSeconds * timeUnitInSeconds * timeUnitInSeconds;
		HydrogenMass = HydrogenMassCGS / massUnitInGrams;
		GasConstant = BoltzmannConstant / HydrogenMass;
	}
	
};

#endif //CODEUNITS_H__hsdf07j253v89724tvuhi23r2t3vuj23v
