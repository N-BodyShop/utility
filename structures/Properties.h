/** \file Properties.h
 * The Properties class represents a set of key/value pairs of different types.
 * It can be read from and written to streams.
 * \author Graeme Lufkin gwl@u.washington.edu
 */
 
#ifndef PROPERTIES_H
#define PROPERTIES_H

#include <fstream>
#include <map>
#include <string>
#include <sstream>
#include <algorithm>

using std::string;
using std::map;
using std::istringstream;
	
class Properties {
private:
	map<string, bool> boolProperties;
	map<string, int> intProperties;
	map<string, double> doubleProperties;
	map<string, string> stringProperties;
	
	/// Remove any leading and trailing whitespace from a string
	static string trim(const string& s) {
		string trimmed(s);
		trimmed.erase(0, trimmed.find_first_not_of(" \t\r\n"));
		trimmed.erase(trimmed.find_last_not_of(" \t\r\n") + 1);
		if(trimmed[trimmed.length() - 1] == ';') {
			trimmed.erase(trimmed.length() - 1);
			trimmed.erase(trimmed.find_last_not_of(" \t\r\n") + 1);
		}
		return trimmed;
	}
	
	/// Remove quotes from around a string, if they are the first and last characters
	static string removeQuotes(const string& s) {
		if((s[0] == '\"' && s[s.length() - 1] == '\"') || (s[0] == '\'' && s[s.length() - 1] == '\''))
			return s.substr(1, s.length() - 2);
		else
			return s;
	}
	
	/// Write this set of key/value pairs to stream
	void store(std::ostream& os) const {
		os << "#bool properties\n";
		for(map<string, bool>::const_iterator iter = boolProperties.begin(), end = boolProperties.end(); iter != end; ++iter)
			if(iter->second)
				os << iter->first << "=true\n";
			else
				os << iter->first << "=false\n";
		os << "#integer properties\n";
		for(map<string, int>::const_iterator intIter = intProperties.begin(); intIter != intProperties.end(); ++intIter)
			os << intIter->first << "=" << intIter->second << "\n";
		os << "#floating point properties\n";
		os.setf(std::ios_base::scientific);
		for(map<string, double>::const_iterator doubleIter = doubleProperties.begin(); doubleIter != doubleProperties.end(); ++doubleIter)
			os << doubleIter->first << "=" << doubleIter->second << "\n";
		os << "#string properties\n";
		for(map<string, string>::const_iterator iter = stringProperties.begin(), end = stringProperties.end(); iter != end; ++iter)
			os << iter->first << "=\"" << iter->second << "\"\n";
	}
	
public:
	
	/// Create an empty Properties
	Properties() { }
	
	/// Create a Properties from an input stream
	Properties(std::istream& in) {
		load(in);
	}
	
	/// Create a Properties from a file
	Properties(const char* filename) {
		std::ifstream in(filename);
		load(in);
		in.close();
	}
	
	/// Copy a Properties
	Properties(const Properties& p) {
		boolProperties = p.boolProperties;
		intProperties = p.intProperties;
		doubleProperties = p.doubleProperties;
		stringProperties = p.stringProperties;
	}
	
	~Properties() { }
	
	Properties& operator=(const Properties& p) {
		boolProperties = p.boolProperties;
		intProperties = p.intProperties;
		doubleProperties = p.doubleProperties;
		stringProperties = p.stringProperties;
		return *this;
	}
	
	/// Two Properties are equal if all their key/value pairs match
	bool operator==(const Properties& p) const {
		return boolProperties == p.boolProperties
				&& intProperties == p.intProperties 
				&& doubleProperties == p.doubleProperties 
				&& stringProperties == p.stringProperties;
	}
	
	bool operator!=(const Properties& p) const {
		return boolProperties != p.boolProperties
				|| intProperties != p.intProperties 
				|| doubleProperties != p.doubleProperties 
				|| stringProperties != p.stringProperties;
	}
	
	/// Add another set of Properties to this one, overriding any already existing key/value pairs
	Properties& operator+=(const Properties& p) {
		map<string, bool> newBoolProperties(p.boolProperties);
		newBoolProperties.insert(boolProperties.begin(), boolProperties.end());
		swap(boolProperties, newBoolProperties);
		map<string, int> newIntProperties(p.intProperties);
		newIntProperties.insert(intProperties.begin(), intProperties.end());
		swap(intProperties, newIntProperties);
		map<string, double> newDoubleProperties(p.doubleProperties);
		newDoubleProperties.insert(doubleProperties.begin(), doubleProperties.end());
		swap(doubleProperties, newDoubleProperties);
		map<string, string> newStringProperties(p.stringProperties);
		newStringProperties.insert(stringProperties.begin(), stringProperties.end());
		swap(stringProperties, newStringProperties);
		return *this;
	}
	
	/// Erase all the key/value pairs
	void clear() {
		boolProperties.clear();
		intProperties.clear();
		doubleProperties.clear();
		stringProperties.clear();
	}
	
	bool getBoolProperty(const string& key, const bool fallback = false) const {
		map<string, bool>::const_iterator iter = boolProperties.find(key);
		if(iter != boolProperties.end())
			return iter->second;
		else
			return fallback;
	}
	
	int getIntProperty(const string& key, const int fallback = 0) const {
		map<string, int>::const_iterator iter = intProperties.find(key);
		if(iter != intProperties.end())
			return iter->second;
		else
			return fallback;
	}
	
	double getDoubleProperty(const string& key, const double fallback = 0.0) const {
		map<string, double>::const_iterator iter = doubleProperties.find(key);
		if(iter != doubleProperties.end())
			return iter->second;
		else
			return fallback;
	}
	
	string getStringProperty(const string& key, const string fallback = "") const {
		map<string, string>::const_iterator iter = stringProperties.find(key);
		if(iter != stringProperties.end())
			return iter->second;
		else
			return fallback;
	}
	
	void setProperty(const string& key, const bool val) {
		boolProperties[key] = val;
	}
	
	void setProperty(const string& key, const int val) {
		intProperties[key] = val;
	}

	void setProperty(const string& key, const double val) {
		doubleProperties[key] = val;
	}

	void setProperty(const string& key, const string val) {
		stringProperties[key] = val;
	}
	
	/// Needed because (apparently!) const char* to bool is a closer match than const char* to string
	void setProperty(const string& key, const char* val) {
		stringProperties[key] = val;
	}
	
	bool eraseBoolProperty(const string& key) {
		return boolProperties.erase(key);
	}
	
	bool eraseIntProperty(const string& key) {
		return intProperties.erase(key);
	}
	
	bool eraseDoubleProperty(const string& key) {
		return intProperties.erase(key);
	}
	
	bool eraseStringProperty(const string& key) {
		return stringProperties.erase(key);
	}
	
	/// Remove a key and its value(s) from any of the maps
	bool eraseProperty(const string& key) {
		return eraseBoolProperty(key)
				|| eraseIntProperty(key)
				|| eraseDoubleProperty(key)
				|| eraseStringProperty(key);
	}
	
	void load(const std::string& s) {
		std::ifstream infile(s.c_str());
		if(infile)
			load(infile);
		infile.close();	
	}
	
	/// Read in a set of key/value pairs from a stream, recognizing types and skipping '#' comments
	void load(std::istream& in) {
		string line, key, value;
		std::size_t equalsPos, commentPos;
		double doubleVal;
		int intVal;
		for(;;) {
			getline(in, line); //read in a line
			if(!in) //check if we're done
				return;
			equalsPos = line.find('='); //look for equals sign
			if(equalsPos > 0) { //if there's an equals
				commentPos = line.find('#'); //look for a comment start
				if(commentPos == string::npos) //pretend comment at end if none
					commentPos = line.size();
				if(commentPos > equalsPos) { //make sure the equals sign is in the right place
					key = trim(line.substr(0, equalsPos));
					value = trim(line.substr(equalsPos + 1, commentPos - equalsPos - 1));
					istringstream intIn(value);
					intIn >> intVal;
					if(intIn && (intIn.get() == -1))
						intProperties[key] = intVal;
					else {
						istringstream doubleIn(value);
						doubleIn >> doubleVal;
						if(doubleIn && (doubleIn.get() == -1)) //successful extraction of double
							doubleProperties[key] = doubleVal;
						else { //it's not an int or double
							//look for boolean values
							string lowervalue(value);
							std::transform(value.begin(), value.end(), lowervalue.begin(), tolower);
							if(lowervalue == "true" || lowervalue == "yes" || lowervalue == "on")
								boolProperties[key] = true;
							else if(lowervalue == "false" || lowervalue == "no" || lowervalue == "off")
								boolProperties[key] = false;
							else //it's got to be a string
								stringProperties[key] = removeQuotes(value);
						}
					}
				}
			}
		}
	}
	
	/// Print this set of key/value pairs to stream
	friend std::ostream& operator<<(std::ostream& os, const Properties& p) {
		p.store(os);
		return os;
	}
	
};


#endif //PROPERTIES_H
