/** @file tipsy2tree.cpp
 Convert a tipsy format file of particle data into tree-based
 set of files.
 @author Graeme Lufkin (gwl@u.washington.edu)
 @date Created February 12, 2003
 @version 1.0
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

#include "TipsyReader.h"
#include "OrientedBox.h"
#include "tree_xdr.h"
#include "SFC.h"

using namespace std;
using namespace SFC;
using namespace Tipsy;

/// A particle that contains all the attributes for dark matter particles, plus some
class SFCParticle {
public:
	
	Key key;
	unsigned int tipsyOrder; //the original order on disk
	float mass;
	Vector3D<float> position;
	Vector3D<float> velocity;
	float softening;
	float potential;
	
	SFCParticle(Key k = 0) : key(k) { }

	//All particles get turned into this type, ignoring extra fields
	SFCParticle(const simple_particle& p) : mass(p.mass), position(p.pos), velocity(p.vel), softening(0), potential(0) { }	
	SFCParticle(const gas_particle& p) :  mass(p.mass), position(p.pos), velocity(p.vel), softening(p.hsmooth), potential(p.phi) { }
	SFCParticle(const dark_particle& p) :  mass(p.mass), position(p.pos), velocity(p.vel), softening(p.eps), potential(p.phi) { }
	SFCParticle(const star_particle& p) :  mass(p.mass), position(p.pos), velocity(p.vel), softening(p.eps), potential(p.phi) { }

	
	/** Comparison operator, used for a sort based on key values. */
	bool operator<(const SFCParticle& p) const {
		return key < p.key;
	}
};

enum NodeType {
	Invalid,
	Bucket,
	Internal,
	NonLocal,
	Empty,
	Boundary,
	Top
};

int maxBucketSize;

class MyTreeNode  {
	NodeType myType;
public:
	
	Key key;
	unsigned char level;
	
	unsigned int numNodesLeft;
	unsigned int numNodesRight;
	unsigned int numParticlesLeft;
	unsigned int numParticlesRight;
	
	union {
		MyTreeNode* leftChild;
		SFCParticle* beginBucket;
	};

	union {
		MyTreeNode* rightChild;
		SFCParticle* endBucket;
	};
	
	MyTreeNode() : myType(Invalid), key(0), level(0), numNodesLeft(0), numNodesRight(0), numParticlesLeft(0), numParticlesRight(0), leftChild(0), rightChild(0) { }
	
	MyTreeNode* createLeftChild() {
		MyTreeNode* child = new MyTreeNode;
		child->key = key;
		child->level = level + 1;
		leftChild = child;
		return child;
	}
	
	MyTreeNode* createRightChild() {
		MyTreeNode* child = new MyTreeNode;
		child->key = key | (static_cast<Key>(1) << (62 - level));
		child->level = level + 1;
		rightChild = child;
		return child;
	}
	
	~MyTreeNode() {
		if(myType != Bucket) {
			delete leftChild;
			delete rightChild;
		}
	}
	
	NodeType getType() const { return myType; }
	void setType(NodeType t) { myType = t; }
	
};

SFCParticle* leftBoundary;
SFCParticle* rightBoundary;

unsigned int totalNumNodes;
unsigned int totalNumParticles;

unsigned int* bucketCounts;

void buildTree(MyTreeNode* node, SFCParticle* leftParticle, SFCParticle* rightParticle) {
	//check if we should bucket these particles
	if(rightParticle - leftParticle < maxBucketSize) {
		//can't bucket until we've cut at the boundary
		if((leftParticle != leftBoundary) && (rightParticle != rightBoundary)) {
			node->setType(Bucket);
			node->beginBucket = leftParticle;
			node->endBucket = rightParticle + 1;
			totalNumParticles += node->endBucket - node->beginBucket;
			node->numParticlesLeft = node->endBucket - node->beginBucket;
			bucketCounts[node->endBucket - node->beginBucket]++;
			node->numParticlesRight = 0;
			node->numNodesLeft = 0;
			node->numNodesRight = 0;
			return;
		}
	} else if(node->level == 63) {
		cerr << "This tree has exhausted all the bits in the keys.  Super double-plus ungood!" << endl;
		return;
	}
	
	//this is the bit we are looking at
	Key currentBitMask = static_cast<Key>(1) << (62 - node->level);
	//we need to know the bit values at the left and right
	Key leftBit = leftParticle->key & currentBitMask;
	Key rightBit = rightParticle->key & currentBitMask;
	MyTreeNode* child;
	
	node->setType(Internal);
	
	if(leftBit ^ rightBit) { //a split at this level
		//find the split by looking for where the key with the bit switched on could go
		SFCParticle* splitParticle = lower_bound(leftParticle, rightParticle + 1, node->key | currentBitMask);
		if(splitParticle == leftBoundary + 1) {
			bucketCounts[0]++;
			child = node->createRightChild();
			buildTree(child, splitParticle, rightParticle);
		} else if(splitParticle == rightBoundary) {
			child = node->createLeftChild();
			buildTree(child, leftParticle, splitParticle - 1);
			bucketCounts[0]++;
		} else {
			child = node->createLeftChild();
			buildTree(child, leftParticle, splitParticle - 1);
			child = node->createRightChild();
			buildTree(child, splitParticle, rightParticle);
		}
	} else if(leftBit & rightBit) { //both ones, make a right child
		bucketCounts[0]++;
		child = node->createRightChild();
		buildTree(child, leftParticle, rightParticle);
	} else { //both zeros, make a left child
		child = node->createLeftChild();
		buildTree(child, leftParticle, rightParticle);
		bucketCounts[0]++;
	}
	
	totalNumNodes++;
	if(node->leftChild) {
		node->numParticlesLeft = node->leftChild->numParticlesLeft + node->leftChild->numParticlesRight;
		node->numNodesLeft = node->leftChild->numNodesLeft + node->leftChild->numNodesRight + (node->leftChild->getType() == Bucket ? 0 : 1);
	} else {
		node->numParticlesLeft = 0;
		node->numNodesLeft = 0;
	}
	if(node->rightChild) {
		node->numParticlesRight = node->rightChild->numParticlesLeft + node->rightChild->numParticlesRight;
		node->numNodesRight = node->rightChild->numNodesLeft + node->rightChild->numNodesRight + (node->rightChild->getType() == Bucket ? 0 : 1);
	} else {
		node->numParticlesRight = 0;
		node->numNodesRight = 0;
	}
}

int xdr_convert_tree(XDR* xdrs, MyTreeNode* node) {
	if(node->getType() != Bucket) {
		u_int64_t numNodesLeft = node->numNodesLeft;
		u_int64_t numParticlesLeft = node->numParticlesLeft;
		int result = xdr_template(xdrs, &numNodesLeft) && xdr_template(xdrs, &numParticlesLeft);
		if(node->leftChild)
			result = result && xdr_convert_tree(xdrs, node->leftChild);
		if(node->rightChild)
			result = result && xdr_convert_tree(xdrs, node->rightChild);
		return result;
	} else
		return 1;
}

void print_tree(OrientedBox<double> box, MyTreeNode* node, int axis) {
	if(node->getType() == Bucket)
		return;
	cout << "Box: " << box << endl;
	cout << "Node: numNodesLeft: " << node->numNodesLeft << " numNodesRight: " << node->numNodesRight << endl;
	cout << "numParticlesLeft: " << node->numParticlesLeft << " numParticlesRight: " << node->numParticlesRight << endl;
	if(node->leftChild)
		print_tree(cutBoxLeft(box, axis), node->leftChild, (axis + 1) % 3);
	if(node->rightChild)
		print_tree(cutBoxRight(box, axis), node->rightChild, (axis + 1) % 3);
}

int main(int argc, char** argv) {

	if(argc < 2) {
		cerr << "Usage: " << argv[0] << " tipsyfile [bucketSize]" << endl;
		return 1;
	}
	
	cerr << "Loading tipsy file ..." << endl;
	TipsyReader r(argv[1]);
	if(!r.status()) {
		cerr << "Couldn't open tipsy file correctly!  Maybe it's not a tipsy file?" << endl;
		return 2;
	}
	
	if(argc > 2) {
		maxBucketSize = atoi(argv[2]);
		if(maxBucketSize <= 0) {
			maxBucketSize = 12;
			cerr << "Bucket size must be positive, using default value (" << maxBucketSize << ")" << endl;
		}
	} else
		maxBucketSize = 12;
	
	header h = r.getHeader();
	SFCParticle* particles = new SFCParticle[h.nbodies + 2];
	gas_particle gp;
	dark_particle dp;
	star_particle sp;
	cerr << "Tipsy header:\n" << h << endl;
	
	OrientedBox<Real> boundingBox;
	cerr << "Converting mass, position, velocity, softening, and potential" << endl;
	for(unsigned int i = 0; i < h.nsph; ++i) {
		if(!r.getNextGasParticle(gp)) {
			cerr << "Error reading tipsy file!" << endl;
			return 3;
		}
		boundingBox.grow(gp.pos);
		particles[i] = gp;
		particles[i].tipsyOrder = i;
	}
	for(unsigned int i = 0; i < h.ndark; ++i) {
		if(!r.getNextDarkParticle(dp)) {
			cerr << "Error reading tipsy file!" << endl;
			return 3;
		}
		boundingBox.grow(dp.pos);
		particles[i] = dp;
		particles[i].tipsyOrder = i + h.nsph;
	}
	for(unsigned int i = 0; i < h.nstar; ++i) {
		if(!r.getNextStarParticle(sp)) {
			cerr << "Error reading tipsy file!" << endl;
			return 3;
		}
		boundingBox.grow(sp.pos);
		particles[i] = sp;
		particles[i].tipsyOrder = i + h.nsph + h.ndark;
	}
	
	OrientedBox<Real> unitCube(Vector3D<Real>(-0.5, -0.5, -0.5), Vector3D<Real>(0.5, 0.5, 0.5));
	Real diagonal = (boundingBox.greater_corner - boundingBox.lesser_corner).length();
	Real epsilon = max(1E-2, 1.0 / h.nbodies);
	if((boundingBox.lesser_corner - unitCube.lesser_corner).length() / diagonal < epsilon
			&& (boundingBox.greater_corner - unitCube.greater_corner).length() / diagonal < epsilon) {
		cerr << "The bounding box for this file appears to be a unit cube about the origin, using this for key generation" << endl;
		boundingBox = unitCube;
	} else {
		cerr << "The bounding box for this file is " << boundingBox << ", using longest axis to generate keys" << endl;
		cubize(boundingBox);
		bumpBox(boundingBox, HUGE_VAL);
		cerr << "Resized bounding box is " << boundingBox << endl;
	}
	
	cerr << "Calculating keys ..." << endl;
	for(unsigned int i = 0; i < h.nbodies; ++i)
		particles[i].key = generateKey(particles[i].position, boundingBox);
	
	cerr << "Sorting particles ..." << endl;
	SFCParticle dummy;
	dummy.key = firstPossibleKey;
	particles[h.nbodies] = dummy;
	dummy.key = lastPossibleKey;
	particles[h.nbodies + 1] = dummy;
	sort(particles, particles + h.nbodies + 2);
	leftBoundary = particles;
	rightBoundary = particles + h.nbodies + 2 - 1;
		
	//build tree
	bucketCounts = new unsigned int[maxBucketSize + 1];
	for(int i = 0; i <= maxBucketSize; ++i)
		bucketCounts[i] = 0;
	cerr << "Building the tree ..." << endl;
	MyTreeNode* root = new MyTreeNode;
	root->key = firstPossibleKey;
	buildTree(root, leftBoundary, rightBoundary);
	
	cerr << "Total number of internal nodes: " << totalNumNodes << endl;
	cerr << "Total number of particles: " << totalNumParticles << endl;
	unsigned int numBuckets = 0;
	for(int i = 1; i <= maxBucketSize; ++i)
		numBuckets += bucketCounts[i];
	cerr << "Total number of leaf nodes: " << (numBuckets + bucketCounts[0]) << endl;
	cerr << "Total number of buckets: " << numBuckets << endl;
	cerr << "Bucket statistics:" << endl;
	for(int i = 1; i <= maxBucketSize; ++i)
		cerr << i << "\t" << ((double) bucketCounts[i] / numBuckets) << endl;
	cerr << "Average particles per bucket: " << ((double) totalNumParticles / numBuckets) << endl;
	delete[] bucketCounts;
	
	TreeHeader th;
	th.time = h.time;
	th.numNodes = totalNumNodes;
	th.numParticles = totalNumParticles;
	th.boundingBox = boundingBox;
	
	string basename(argv[1]);
	int slashPos = basename.find_last_of("/");
	if(slashPos != string::npos)
		basename.erase(0, slashPos + 1);
	cerr << "Basename for tree-based files: \"" << basename << "\"" << endl;
	
	//write xml
	ofstream xmlfile((basename + ".xml").c_str());
	xmlfile << "<xml fakery>\n<simulation>\n";
	xmlfile << "\t<tree filename=\"" << basename << ".tree\">\n";
	xmlfile << "\t<mass filename=\"" << basename << ".mass\">\n";
	xmlfile << "\t<position filename=\"" << basename << ".pos\">\n";
	xmlfile << "\t<velocity filename=\"" << basename << ".vel\">\n";
	xmlfile << "</simulation>\n";
	xmlfile.close();
	
	FILE* outfile;
	XDR xdrs;

	//write tree file
	outfile = fopen((basename + ".tree").c_str(), "wb");
	xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
	xdr_template(&xdrs, &th);
	xdr_convert_tree(&xdrs, root);
	xdr_destroy(&xdrs);
	fclose(outfile);
	
	//print_tree(th.boundingBox, root, 0);
	
	delete root;
	
	FieldHeader fh;
	fh.time = th.time;
	fh.numParticles = totalNumParticles;
	
	//write extras files
	
	outfile = fopen((basename + ".uid").c_str(), "wb");
	xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
	fh.dimensions = 1;
	fh.code = uint32;
	xdr_template(&xdrs, &fh);
	for(unsigned int i = 1; i <= totalNumParticles; ++i)
		xdr_template(&xdrs, &(particles[i].tipsyOrder));
	xdr_destroy(&xdrs);
	fclose(outfile);

	outfile = fopen((basename + ".mass").c_str(), "wb");
	xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
	fh.dimensions = 1;
	fh.code = float32;
	xdr_template(&xdrs, &fh);
	for(unsigned int i = 1; i <= totalNumParticles; ++i)
		xdr_template(&xdrs, &(particles[i].mass));
	xdr_destroy(&xdrs);
	fclose(outfile);
	
	outfile = fopen((basename + ".pos").c_str(), "wb");
	xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
	fh.dimensions = 3;
	fh.code = float32;
	xdr_template(&xdrs, &fh);
	for(unsigned int i = 1; i <= totalNumParticles; ++i)
		xdr_template(&xdrs, &(particles[i].position));
	xdr_destroy(&xdrs);
	fclose(outfile);
	
	outfile = fopen((basename + ".vel").c_str(), "wb");
	xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
	fh.dimensions = 3;
	fh.code = float32;
	xdr_template(&xdrs, &fh);
	for(unsigned int i = 1; i <= totalNumParticles; ++i)
		xdr_template(&xdrs, &(particles[i].velocity));
	xdr_destroy(&xdrs);
	fclose(outfile);
	
	outfile = fopen((basename + ".softening").c_str(), "wb");
	xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
	fh.dimensions = 1;
	fh.code = float32;
	xdr_template(&xdrs, &fh);
	for(unsigned int i = 1; i <= totalNumParticles; ++i)
		xdr_template(&xdrs, &(particles[i].softening));
	xdr_destroy(&xdrs);
	fclose(outfile);
	
	outfile = fopen((basename + ".potential").c_str(), "wb");
	xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
	fh.dimensions = 1;
	fh.code = float32;
	xdr_template(&xdrs, &fh);
	for(unsigned int i = 1; i <= totalNumParticles; ++i)
		xdr_template(&xdrs, &(particles[i].potential));
	xdr_destroy(&xdrs);
	fclose(outfile);

	cerr << "Tipsy format file \"" << basename << "\" converted to tree-based format." << endl;
	
	delete[] particles;
}
