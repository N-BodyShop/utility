//Properties.h

#ifndef PROPERTIES_H
#define PROPERTIES_H

#include <iostream>
#include <map>
#include <string>
#include <sstream>

using std::map;
using std::string;
using std::istringstream;

class Properties {
private:
	map<string, int> intProperties;
	map<string, double> doubleProperties;
	map<string, string> stringProperties;
	
	//the maximum length of a line in a stored Properties
	static const int lineLength = 1024;
	
	void trim(string& s) {
		s.erase(0, s.find_first_not_of(" \t\r\n"));
		s.erase(s.find_last_not_of(" \t\r\n") + 1);
	}
	
public:
	Properties() { }
	
	//create a Properties from an input stream
	Properties(std::istream& in) {
		load(in);
	}
	
	//copy constructor
	Properties(const Properties& p) {
		intProperties = p.intProperties;
		doubleProperties = p.doubleProperties;
		stringProperties = p.stringProperties;
	}
	
	~Properties() { }
	
	void clear() {
		intProperties.clear();
		doubleProperties.clear();
		stringProperties.clear();
	}
	
	bool operator==(const Properties& p) {
		return intProperties == p.intProperties 
				&& doubleProperties == p.doubleProperties 
				&& stringProperties == p.stringProperties;
	}
	
	int getPropertyInt(const string prop) const {
		map<string, int>::const_iterator iter = intProperties.find(prop);
		if(iter != intProperties.end())
			return iter->second;
		else
			return 0;
	}
	
	double getPropertyDouble(const string prop) const {
		map<string, double>::const_iterator iter = doubleProperties.find(prop);
		if(iter != doubleProperties.end())
			return iter->second;
		else
			return 0.0;
	}
	
	string getPropertyString(const string prop) const {
		map<string, string>::const_iterator iter = stringProperties.find(prop);
		if(iter != stringProperties.end())
			return iter->second;
		else
			return "";
	}
	
	void setProperty(const string prop, const int val) {
		intProperties[prop] = val;
	}

	void setProperty(const string prop, const double val) {
		doubleProperties[prop] = val;
	}

	void setProperty(const string prop, const string val) {
		stringProperties[prop] = val;
	}

	void load(std::istream& in) {
		char line[lineLength];
		string s, key, value;
		int equalsPos, poundPos;
		double doubleVal;
		int intVal;
		istringstream iss;
		for(;;) {
			in.getline(line, lineLength); //read in a line
			if(!in) //check if we're done
				return;
			s = line;
			equalsPos = s.find('='); //look for equals sign
			if(equalsPos > 0) { //if there's an equals
				poundPos = s.find('#'); //look for a comment start
				if(poundPos == string::npos) //if there's no pound sign, pretend it's at the end
					poundPos = s.size();
				if(poundPos > equalsPos) { //make sure the equals sign is in the right place
					key = s.substr(0, equalsPos);
					trim(key);
					value = s.substr(equalsPos + 1, poundPos - equalsPos - 1);
					trim(value);
					iss.str(value);
					iss.seekg(0); //make sure we're at the beginning of the stream
					iss.clear(); //clear any error bits
					iss >> intVal;
					if(iss.eof()) //it was an int
						intProperties[key] = intVal;
					else { //try a double
						iss.seekg(0);
						iss.clear();
						iss >> doubleVal;
						if(iss.eof()) //successful extraction of double
							doubleProperties[key] = doubleVal;
						else //it's not an int or double, must be string
							stringProperties[key] = value;
					}
				}			
			}
		}
	}
	
	void store(std::ostream& os) const {
		for(map<string, int>::const_iterator intIter = intProperties.begin(); intIter != intProperties.end(); ++intIter)
			os << intIter->first << "=" << intIter->second << "\n";
		for(map<string, double>::const_iterator doubleIter = doubleProperties.begin(); doubleIter != doubleProperties.end(); ++doubleIter)
			os << doubleIter->first << "=" << doubleIter->second << "\n";
		for(map<string, string>::const_iterator stringIter = stringProperties.begin(); stringIter != stringProperties.end(); ++stringIter)
			os << stringIter->first << "=" << stringIter->second << "\n";
	}
	
	friend std::ostream& operator<< (std::ostream& os, const Properties& p) {
		p.store(os);
		return os;
	}
	
};

#endif //PROPERTIES_H
