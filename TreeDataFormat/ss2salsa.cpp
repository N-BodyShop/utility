/** @file ss2salsa.cpp
 Convert a ss format file of particle data.
 @author Graeme Lufkin (gwl@u.washington.edu)
 @date Created February 13, 2004
 @version 1.0
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <popt.h>

#include "SS.h"
#include "OrientedBox.h"
#include "tree_xdr.h"
#include "aggregate_xdr.h"
#include "SFC.h"

using namespace std;
using namespace SS;
using namespace TypeHandling;
using namespace SFC;

int verbosity;

MAKE_AGGREGATE_WRITER(mass)
MAKE_AGGREGATE_WRITER(radius)
//MAKE_AGGREGATE_WRITER(pos)
MAKE_AGGREGATE_WRITER(vel)
MAKE_AGGREGATE_WRITER(spin)
MAKE_AGGREGATE_WRITER(color)
MAKE_AGGREGATE_WRITER(org_idx)

template <typename T>
void determineBoundingBox(OrientedBox<T>& boundingBox, const int numParticles) {
	cubize(boundingBox);
	//bumpBox(boundingBox, HUGE_VAL);
	if(verbosity > 1)
		cerr << "Resized bounding box is " << boundingBox << endl;
}

bool convertParticles(const string& filenamePrefix, SSReader& r) {
	int numParticles = r.getHeader().n_data;
	
	vector<Particle> particles;
	particles.reserve(numParticles);
	
	r.readAllParticles(back_inserter(particles));
	
	SSStats stats;
	for(int i = 0; i < numParticles; ++i)
		stats.contribute(particles[i]);
	stats.finalize();
	OrientedBox<float> boundingBox(stats.bounding_box);

	determineBoundingBox(boundingBox, numParticles);
	
	//make ss subdirectory
	if(mkdir("ss", 0755))
		return false;
	
	//write out uid, mass, radius, pos, vel, spin, color
	FILE* outfile;
	XDR xdrs;
	
	FieldHeader fh;
	fh.time = r.getHeader().time;
	fh.numParticles = numParticles;
	
	fh.dimensions = 1;
	fh.code = uint32;
	outfile = fopen("ss/uid", "wb");
	xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
	writeAggregateMember_org_idx(&xdrs, fh, &(*particles.begin()), stats.min_org_idx, stats.max_org_idx);
	xdr_destroy(&xdrs);
	fclose(outfile);	
	
	fh.dimensions = 1;
	fh.code = float64;
	outfile = fopen("ss/Mass", "wb");
	xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
	writeAggregateMember_mass(&xdrs, fh, &(*particles.begin()), stats.min_mass, stats.max_mass);
	xdr_destroy(&xdrs);
	fclose(outfile);	

	fh.dimensions = 1;
	fh.code = float64;
	outfile = fopen("ss/Radius", "wb");
	xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
	writeAggregateMember_radius(&xdrs, fh, &(*particles.begin()), stats.min_radius, stats.max_radius);
	xdr_destroy(&xdrs);
	fclose(outfile);	
	
	//currently salsa needs position to be float, not double!
	vector<Vector3D<float> > positions;
	positions.reserve(numParticles);
	for(int i = 0; i < numParticles; ++i)
		positions.push_back(particles[i].pos);
	fh.dimensions = 3;
	fh.code = float32;
	outfile = fopen("ss/position", "wb");
	xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
	writeField(fh, &xdrs, &(*positions.begin()), &boundingBox.lesser_corner, &boundingBox.greater_corner);
	xdr_destroy(&xdrs);
	fclose(outfile);	

	fh.dimensions = 3;
	fh.code = float64;
	outfile = fopen("ss/Velocity", "wb");
	xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
	writeAggregateMember_vel(&xdrs, fh, &(*particles.begin()), stats.velocityBox.lesser_corner, stats.velocityBox.greater_corner);
	xdr_destroy(&xdrs);
	fclose(outfile);	

	fh.dimensions = 3;
	fh.code = float64;
	outfile = fopen("ss/Spin", "wb");
	xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
	writeAggregateMember_vel(&xdrs, fh, &(*particles.begin()), stats.spinBox.lesser_corner, stats.spinBox.greater_corner);
	xdr_destroy(&xdrs);
	fclose(outfile);	

	fh.dimensions = 1;
	fh.code = uint32;
	outfile = fopen("ss/Color", "wb");
	xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
	writeAggregateMember_color(&xdrs, fh, &(*particles.begin()), stats.min_color, stats.max_color);
	xdr_destroy(&xdrs);
	fclose(outfile);	
	
	return true;
}

int main(int argc, const char** argv) {

	poptOption optionsTable[] = {
		{"verbose", 'v', POPT_ARG_NONE | POPT_ARGFLAG_ONEDASH | POPT_ARGFLAG_SHOW_DEFAULT, 0, 1, "be verbose about what's going on", "verbosity"},
		POPT_AUTOHELP
		POPT_TABLEEND
	};
	
	poptContext context = poptGetContext("ss2salsa", argc, argv, optionsTable, 0);
	
	poptSetOtherOptionHelp(context, " [OPTION ...] ssfile");
	
	int rc;
	while((rc = poptGetNextOpt(context)) >= 0) {
		switch(rc) {
			case 1: //increase verbosity
				verbosity++;
				break;
		}
	}
	
	if(rc < -1) {
		cerr << "Argument error: " << poptBadOption(context, POPT_BADOPTION_NOALIAS) << " : " << poptStrerror(rc) << endl;
		poptPrintUsage(context, stderr, 0);
		return 1;
	}
	
	const char* fname = poptGetArg(context);
	
	if(fname == 0) {
		cerr << "You must provide a SS file name" << endl;
		poptPrintUsage(context, stderr, 0);
		return 2;
	}
	
	if(verbosity > 1)
		cerr << "Loading SS file ..." << endl;
	
	SSReader r(fname);
	if(!r.status()) {
		cerr << "Couldn't open SS file correctly!  Maybe it's not a SS file?" << endl;
		return 2;
	}
		
	string basename(fname);
	string::size_type slashPos = basename.find_last_of("/");
	if(slashPos != string::npos)
		basename.erase(0, slashPos + 1);
	cerr << "Basename for tree-based files: \"" << basename << "\"" << endl;
	string dirname = basename + ".data";
	
	//make directory for files
	if(mkdir(dirname.c_str(), 0775) || chdir(dirname.c_str())) {
		cerr << "Could not create and move into a directory for the converted files, maybe you don't have permission?" << endl;
		return 3;
	}
	//make output a string for output index
	
	//write xml
	ofstream xmlfile("description.xml");
	xmlfile << "<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n<simulation>\n";

	Header h = r.getHeader();
	if(verbosity > 2)
		cout << "SS header:\n" << h << endl;
	
	convertParticles(basename, r);
	xmlfile << "\t<family name=\"ss\">\n";
	xmlfile << "\t\t<attribute name=\"uid\" link=\"ss/uid\"/>\n";
	xmlfile << "\t\t<attribute name=\"Mass\" link=\"ss/Mass\"/>\n";
	xmlfile << "\t\t<attribute name=\"Radius\" link=\"ss/Radius\"/>\n";
	xmlfile << "\t\t<attribute name=\"position\" link=\"ss/position\"/>\n";
	xmlfile << "\t\t<attribute name=\"Velocity\" link=\"ss/Velocity\"/>\n";
	xmlfile << "\t\t<attribute name=\"Spin\" link=\"ss/Spin\"/>\n";
	xmlfile << "\t\t<attribute name=\"Color\" link=\"ss/Color\"/>\n";
	xmlfile << "\t</family>\n";
		
	xmlfile << "</simulation>\n";
	xmlfile.close();
	
	cerr << "Done." << endl;
}
