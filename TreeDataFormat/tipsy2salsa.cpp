/** @file tipsy2tree.cpp
 Convert a tipsy format file of particle data into the nchilada format
 set of files without putting them into tree order.
 Based on Graeme Lufkin's tipsy2tree.cpp code.
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

#include "TipsyFile.h"
#include "OrientedBox.h"
#include "tree_xdr.h"
#include "aggregate_xdr.h"

using namespace std;
using namespace Tipsy;
using namespace TypeHandling;

int verbosity;

MAKE_AGGREGATE_WRITER(mass)
MAKE_AGGREGATE_WRITER(pos)
MAKE_AGGREGATE_WRITER(vel)
MAKE_AGGREGATE_WRITER(eps)
MAKE_AGGREGATE_WRITER(phi)
MAKE_AGGREGATE_WRITER(rho)
MAKE_AGGREGATE_WRITER(temp)
MAKE_AGGREGATE_WRITER(hsmooth)
MAKE_AGGREGATE_WRITER(metals)
MAKE_AGGREGATE_WRITER(tform)

// Below are two special cases: OxMassFrac and FeMassFrac will be
// written from metals

template <typename Aggregate, typename ValueType>
bool writeAggregateMember_OxMassFrac(XDR* xdrs, FieldHeader& fh, Aggregate* array, ValueType minValue, ValueType maxValue) {
    if(!xdr_template(xdrs, &fh) || !xdr_template(xdrs, &minValue) || !xdr_template(xdrs, &maxValue))
        return false;
    for(u_int64_t i = 0; i < fh.numParticles; ++i) {
        // O and Fe ratio based on Asplund et al 2009
        ValueType dOxMassFrac = array[i].metals*0.43;
        if(!xdr_template(xdrs, &dOxMassFrac))
            return false;
    }
    return true;
}
template <typename Aggregate, typename ValueType>
bool writeAggregateMember_FeMassFrac(XDR* xdrs, FieldHeader& fh, Aggregate* array, ValueType minValue, ValueType maxValue) {
    if(!xdr_template(xdrs, &fh) || !xdr_template(xdrs, &minValue) || !xdr_template(xdrs, &maxValue))
        return false;
    for(u_int64_t i = 0; i < fh.numParticles; ++i) {
        // O and Fe ratio based on Asplund et al 2009
        ValueType dFeMassFrac = array[i].metals*0.098;
        if(!xdr_template(xdrs, &dFeMassFrac))
            return false;
    }
    return true;
}

// For all of the following "0" means use the tipsy header for the
// particle count.
int64_t iNDark = 0;              // Input ndark for oversize tipsy file.
int64_t iNSph = 0;               // Input nSPH for oversize tipsy file.
int64_t iNStar = 0;              // Input nStar for oversize tipsy file.

// The following xdr_template()s are copied from TipsyReader.cpp
template <typename TPos, typename TVel>
inline bool_t xdr_template(XDR* xdrs, Tipsy::gas_particle_t<TPos,TVel>* p) {
	return (xdr_template(xdrs, &(p->mass))
		&& xdr_template(xdrs, &(p->pos))
		&& xdr_template(xdrs, &(p->vel))
		&& xdr_template(xdrs, &(p->rho))
		&& xdr_template(xdrs, &(p->temp))
		&& xdr_template(xdrs, &(p->hsmooth))
		&& xdr_template(xdrs, &(p->metals))
		&& xdr_template(xdrs, &(p->phi)));
}

template <typename TPos, typename TVel>
inline bool_t xdr_template(XDR* xdrs, Tipsy::dark_particle_t<TPos,TVel>* p) {
	return (xdr_template(xdrs, &(p->mass))
		&& xdr_template(xdrs, &(p->pos))
		&& xdr_template(xdrs, &(p->vel))
		&& xdr_template(xdrs, &(p->eps))
		&& xdr_template(xdrs, &(p->phi)));
}

template <typename TPos, typename TVel>
inline bool_t xdr_template(XDR* xdrs, Tipsy::star_particle_t<TPos,TVel>* p) {
	return (xdr_template(xdrs, &(p->mass))
		&& xdr_template(xdrs, &(p->pos))
		&& xdr_template(xdrs, &(p->vel))
		&& xdr_template(xdrs, &(p->metals))
		&& xdr_template(xdrs, &(p->tform))
		&& xdr_template(xdrs, &(p->eps))
		&& xdr_template(xdrs, &(p->phi)));
}
template <typename Part>
void readTipsyParticle(TipsyReader& r, Part *p) 
{
    if(r.isNative()) {
        r.tipsyStream->read((char *)p, sizeof(*p));
    }
    else {
        XDR xdrs;
        char buf[p->sizeBytes];
        r.tipsyStream->read(buf, p->sizeBytes);
        xdrmem_create(&xdrs, buf, p->sizeBytes, XDR_DECODE);
        if(!xdr_template(&xdrs, p))
            assert(0);
        xdr_destroy(&xdrs);
    }
}

bool convertGasParticles(const string& filenamePrefix, TipsyReader& r) {
    int64_t numParticles = r.getHeader().nsph;
    if(iNSph > 0) {
        assert((iNSph & 0xffffffff) == numParticles);
        numParticles = iNSph;
    }

    TipsyStats stats;
    OrientedBox<float> velocityBox;
    OrientedBox<float> boundingBox;

    gas_particle* particles = new gas_particle[numParticles];
    gas_particle p;

    //read in the gas particles from the tipsy file
    for(int64_t i = 0; i < numParticles; ++i) {
        readTipsyParticle(r, &p);
        stats.contribute(p);
        velocityBox.grow(p.vel);
        particles[i] = p;
    }

    stats.finalize();
    boundingBox = stats.bounding_box;

    //make gas subdirectory
    if(mkdir("gas", 0755))
        return false;

    //write out uid, mass, pos, vel, rho, temp, metals, eps, phi
    FILE* outfile;
    XDR xdrs;

    FieldHeader fh;
    fh.time = r.getHeader().time;
    fh.numParticles = numParticles;

    fh.dimensions = 1;
    fh.code = float32;
    outfile = fopen("gas/mass", "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    writeAggregateMember_mass(&xdrs, fh, particles, stats.gas_min_mass, stats.gas_max_mass);
    xdr_destroy(&xdrs);
    fclose(outfile);	

    fh.dimensions = 3;
    fh.code = float32;
    outfile = fopen("gas/pos", "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    writeAggregateMember_pos(&xdrs, fh, particles, boundingBox.lesser_corner, boundingBox.greater_corner);
    xdr_destroy(&xdrs);
    fclose(outfile);	

    fh.dimensions = 3;
    fh.code = float32;
    outfile = fopen("gas/vel", "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    writeAggregateMember_vel(&xdrs, fh, particles, velocityBox.lesser_corner, velocityBox.greater_corner);
    xdr_destroy(&xdrs);
    fclose(outfile);	

    fh.dimensions = 1;
    fh.code = float32;
    outfile = fopen("gas/GasDensity", "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    writeAggregateMember_rho(&xdrs, fh, particles, stats.gas_min_rho, stats.gas_max_rho);
    xdr_destroy(&xdrs);
    fclose(outfile);	

    fh.dimensions = 1;
    fh.code = float32;
    outfile = fopen("gas/temperature", "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    writeAggregateMember_temp(&xdrs, fh, particles, stats.gas_min_temp, stats.gas_max_temp);
    xdr_destroy(&xdrs);
    fclose(outfile);	

    fh.dimensions = 1;
    fh.code = float32;
    outfile = fopen("gas/soft", "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    writeAggregateMember_hsmooth(&xdrs, fh, particles, stats.gas_min_hsmooth, stats.gas_max_hsmooth);
    xdr_destroy(&xdrs);
    fclose(outfile);	

    fh.dimensions = 1;
    fh.code = float32;
    outfile = fopen("gas/metals", "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    writeAggregateMember_metals(&xdrs, fh, particles, stats.gas_min_metals, stats.gas_max_metals);
    xdr_destroy(&xdrs);
    fclose(outfile);	

    fh.dimensions = 1;
    fh.code = float32;
    outfile = fopen("gas/OxMassFrac", "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    writeAggregateMember_OxMassFrac(&xdrs, fh, particles, stats.gas_min_metals, stats.gas_max_metals);
    xdr_destroy(&xdrs);
    fclose(outfile);	

    fh.dimensions = 1;
    fh.code = float32;
    outfile = fopen("gas/FeMassFrac", "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    writeAggregateMember_FeMassFrac(&xdrs, fh, particles, stats.gas_min_metals, stats.gas_max_metals);
    xdr_destroy(&xdrs);
    fclose(outfile);	

    fh.dimensions = 1;
    fh.code = float32;
    outfile = fopen("gas/pot", "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    writeAggregateMember_phi(&xdrs, fh, particles, stats.gas_min_phi, stats.gas_max_phi);
    xdr_destroy(&xdrs);
    fclose(outfile);	

    delete[] particles;

    return true;
}

bool convertDarkParticles(const string& filenamePrefix, TipsyReader& r) {
    int64_t numParticles = r.getHeader().ndark;
    TipsyStats stats;
    OrientedBox<float> velocityBox;
    OrientedBox<float> boundingBox;
	
    if(iNDark > 0) {
        assert((iNDark & 0xffffffff) == numParticles);
        numParticles = iNDark;
    }
    dark_particle* particles = new dark_particle[numParticles];
    dark_particle p;

    //read in the gas particles from the tipsy file
    for(int64_t i = 0; i < numParticles; ++i) {
        readTipsyParticle(r, &p);
        stats.contribute(p);
        velocityBox.grow(p.vel);
        particles[i] = p;
    }

    stats.finalize();
    boundingBox = stats.bounding_box;

    //make dark subdirectory
    if(mkdir("dark", 0775))
            return false;

    //write out uid, mass, pos, vel, eps, phi
    FILE* outfile;
    XDR xdrs;

    FieldHeader fh;
    fh.time = r.getHeader().time;
    fh.numParticles = numParticles;

    fh.dimensions = 1;
    fh.code = float32;
    outfile = fopen("dark/mass", "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    writeAggregateMember_mass(&xdrs, fh, particles, stats.dark_min_mass, stats.dark_max_mass);
    xdr_destroy(&xdrs);
    fclose(outfile);	

    fh.dimensions = 3;
    fh.code = float32;
    outfile = fopen("dark/pos", "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    writeAggregateMember_pos(&xdrs, fh, particles, boundingBox.lesser_corner, boundingBox.greater_corner);
    xdr_destroy(&xdrs);
    fclose(outfile);	

    fh.dimensions = 3;
    fh.code = float32;
    outfile = fopen("dark/vel", "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    writeAggregateMember_vel(&xdrs, fh, particles, velocityBox.lesser_corner, velocityBox.greater_corner);
    xdr_destroy(&xdrs);
    fclose(outfile);	

    fh.dimensions = 1;
    fh.code = float32;
    outfile = fopen("dark/soft", "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    writeAggregateMember_eps(&xdrs, fh, particles, stats.dark_min_eps, stats.dark_max_eps);
    xdr_destroy(&xdrs);
    fclose(outfile);	

    fh.dimensions = 1;
    fh.code = float32;
    outfile = fopen("dark/pot", "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    writeAggregateMember_phi(&xdrs, fh, particles, stats.dark_min_phi, stats.dark_max_phi);
    xdr_destroy(&xdrs);
    fclose(outfile);	

    delete[] particles;

    return true;
}

bool convertStarParticles(const string& filenamePrefix, TipsyReader& r) {
    int64_t numParticles = r.getHeader().nstar;
    TipsyStats stats;
    OrientedBox<float> velocityBox;
    OrientedBox<float> boundingBox;

    if(iNStar > 0) {
        assert((iNStar & 0xffffffff) == numParticles);
        numParticles = iNDark;
    }
    star_particle* particles = new star_particle[numParticles];
    star_particle p;

    //read in the star particles from the tipsy file
    for(int64_t i = 0; i < numParticles; ++i) {
        readTipsyParticle(r, &p);
        stats.contribute(p);
        velocityBox.grow(p.vel);
        particles[i] = p;
    }

    stats.finalize();
    boundingBox = stats.bounding_box;

    //make star subdirectory
    if(mkdir("star", 0775))
        return false;

    //write out uid, mass, pos, vel, metals, tform, eps, phi
    FILE* outfile;
    XDR xdrs;

    FieldHeader fh;
    fh.time = r.getHeader().time;
    fh.numParticles = numParticles;

    fh.dimensions = 1;
    fh.code = float32;
    outfile = fopen("star/mass", "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    writeAggregateMember_mass(&xdrs, fh, particles, stats.star_min_mass, stats.star_max_mass);
    xdr_destroy(&xdrs);
    fclose(outfile);	

    fh.dimensions = 3;
    fh.code = float32;
    outfile = fopen("star/pos", "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    writeAggregateMember_pos(&xdrs, fh, particles, boundingBox.lesser_corner, boundingBox.greater_corner);
    xdr_destroy(&xdrs);
    fclose(outfile);	

    fh.dimensions = 3;
    fh.code = float32;
    outfile = fopen("star/vel", "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    writeAggregateMember_vel(&xdrs, fh, particles, velocityBox.lesser_corner, velocityBox.greater_corner);
    xdr_destroy(&xdrs);
    fclose(outfile);	

    fh.dimensions = 1;
    fh.code = float32;
    outfile = fopen("star/metals", "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    writeAggregateMember_metals(&xdrs, fh, particles, stats.star_min_metals, stats.star_max_metals);
    xdr_destroy(&xdrs);
    fclose(outfile);	

    fh.dimensions = 1;
    fh.code = float32;
    outfile = fopen("star/OxMassFrac", "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    writeAggregateMember_OxMassFrac(&xdrs, fh, particles, stats.gas_min_metals, stats.gas_max_metals);
    xdr_destroy(&xdrs);
    fclose(outfile);	

    fh.dimensions = 1;
    fh.code = float32;
    outfile = fopen("star/FeMassFrac", "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    writeAggregateMember_FeMassFrac(&xdrs, fh, particles, stats.gas_min_metals, stats.gas_max_metals);
    xdr_destroy(&xdrs);
    fclose(outfile);	

    fh.dimensions = 1;
    fh.code = float32;
    outfile = fopen("star/timeform", "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    writeAggregateMember_tform(&xdrs, fh, particles, stats.star_min_tform, stats.star_max_tform);
    xdr_destroy(&xdrs);
    fclose(outfile);	

    fh.dimensions = 1;
    fh.code = float32;
    outfile = fopen("star/soft", "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    writeAggregateMember_eps(&xdrs, fh, particles, stats.star_min_eps, stats.star_max_eps);
    xdr_destroy(&xdrs);
    fclose(outfile);	

    fh.dimensions = 1;
    fh.code = float32;
    outfile = fopen("star/pot", "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    writeAggregateMember_phi(&xdrs, fh, particles, stats.star_min_phi, stats.star_max_phi);
    xdr_destroy(&xdrs);
    fclose(outfile);	

    delete[] particles;

    return true;
}

int main(int argc, char** argv) {

    
    const char *optstring = "d:g:s:v";
    int c;
    while((c=getopt(argc,argv,optstring))>0){
        if(c == -1){
            break;
        }
        switch(c){
        case 'v':
            verbosity++;
            break;
        case 'd':
            iNDark = atoll(optarg);
            break;
        case 'g':
            iNSph = atoll(optarg);
            break;
        case 's':
            iNStar = atoll(optarg);
            break;
        };
    }
    const char *fname;
    if(optind  < argc){
        fname = argv[optind];
    }else{
        cerr << "You must provide a simulation list file" << endl;
        exit(-1);
    }

    if(verbosity > 1)
        cerr << "Loading tipsy file ..." << endl;
	
    TipsyReader r(fname);
    if(!r.status()) {
        cerr << "Couldn't open tipsy file correctly!  Maybe it's not a tipsy file?" << endl;
        return 2;
    }
		
    string basename(fname);
    unsigned int slashPos = basename.find_last_of("/");
    if(slashPos != string::npos)
        basename.erase(0, slashPos + 1);
    cerr << "Basename for tree-based files: \"" << basename << "\"" << endl;
    string dirname = basename + ".data";
	
    //make directory for files
    if(mkdir(dirname.c_str(), 0775) || chdir(dirname.c_str())) {
        cerr << "Could not create and move into a directory for the converted files, maybe you don't have permission?" << endl;
        return 3;
    }

    //write xml
    ofstream xmlfile("description.xml");
    xmlfile << "<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n<simulation>\n";

    header h = r.getHeader();
    if(verbosity > 2)
            cout << "Tipsy header:\n" << h << endl;

    if(h.nsph > 0 || iNSph > 0) {
        convertGasParticles(basename, r);
        xmlfile << "\t<family name=\"gas\">\n";
        xmlfile << "\t\t<attribute name=\"mass\" link=\"gas/mass\"/>\n";
        xmlfile << "\t\t<attribute name=\"position\" link=\"gas/pos\"/>\n";
        xmlfile << "\t\t<attribute name=\"velocity\" link=\"gas/vel\"/>\n";
        xmlfile << "\t\t<attribute name=\"softening\" link=\"gas/soft\"/>\n";
        xmlfile << "\t\t<attribute name=\"potential\" link=\"gas/pot\"/>\n";
        xmlfile << "\t\t<attribute name=\"density\" link=\"gas/GasDensity\"/>\n";
        xmlfile << "\t\t<attribute name=\"temperature\" link=\"gas/temperature\"/>\n";
        xmlfile << "\t\t<attribute name=\"metals\" link=\"gas/metals\"/>\n";
        xmlfile << "\t\t<attribute name=\"OxMassFrac\" link=\"gas/OxMassFrac\"/>\n";
        xmlfile << "\t\t<attribute name=\"FeMassFrac\" link=\"gas/FeMassFrac\"/>\n";
        xmlfile << "\t</family>\n";
    }

    if(h.ndark > 0 || iNDark > 0) {
        convertDarkParticles(basename, r);
        xmlfile << "\t<family name=\"dark\">\n";
        xmlfile << "\t\t<attribute name=\"mass\" link=\"dark/mass\"/>\n";
        xmlfile << "\t\t<attribute name=\"position\" link=\"dark/pos\"/>\n";
        xmlfile << "\t\t<attribute name=\"velocity\" link=\"dark/vel\"/>\n";
        xmlfile << "\t\t<attribute name=\"softening\" link=\"dark/soft\"/>\n";
        xmlfile << "\t\t<attribute name=\"potential\" link=\"dark/pot\"/>\n";
        xmlfile << "\t</family>\n";
    }

    if(h.nstar > 0 || iNStar > 0) {
        convertStarParticles(basename, r);
        xmlfile << "\t<family name=\"star\">\n";
        xmlfile << "\t\t<attribute name=\"mass\" link=\"star/mass\"/>\n";
        xmlfile << "\t\t<attribute name=\"position\" link=\"star/pos\"/>\n";
        xmlfile << "\t\t<attribute name=\"velocity\" link=\"star/vel\"/>\n";
        xmlfile << "\t\t<attribute name=\"softening\" link=\"star/soft\"/>\n";
        xmlfile << "\t\t<attribute name=\"potential\" link=\"star/pot\"/>\n";
        xmlfile << "\t\t<attribute name=\"metals\" link=\"star/metals\"/>\n";
        xmlfile << "\t\t<attribute name=\"OxMassFrac\" link=\"star/OxMassFrac\"/>\n";
        xmlfile << "\t\t<attribute name=\"FeMassFrac\" link=\"star/FeMassFrac\"/>\n";
        xmlfile << "\t\t<attribute name=\"formationtime\" link=\"star/timeform\"/>\n";
        xmlfile << "\t</family>\n";
    }
	
    xmlfile << "</simulation>\n";
    xmlfile.close();

    cerr << "Done." << endl;
}
