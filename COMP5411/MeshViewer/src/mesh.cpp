#include "mesh.h"
#include <iostream>
#include <igl/read_triangle_mesh.h>
#include <Eigen/Sparse>


HEdge::HEdge(bool b) {
	mBoundary = b;

	mTwin = nullptr;
	mPrev = nullptr;
	mNext = nullptr;

	mStart = nullptr;
	mFace = nullptr;

	mFlag = false;
	mValid = true;
}

HEdge* HEdge::twin() const {
	return mTwin;
}

HEdge* HEdge::setTwin(HEdge* e) {
	mTwin = e;
	return mTwin;
}

HEdge* HEdge::prev() const {
	return mPrev;
}

HEdge* HEdge::setPrev(HEdge* e) {
	mPrev = e;
	return mPrev;
}

HEdge* HEdge::next() const {
	return mNext;
}

HEdge* HEdge::setNext(HEdge* e) {
	mNext = e;
	return mNext;
}

Vertex* HEdge::start() const {
	return mStart;
}

Vertex* HEdge::setStart(Vertex* v) {
	mStart = v;
	return mStart;
}

Vertex* HEdge::end() const {
	return mNext->start();
}

Face* HEdge::leftFace() const {
	return mFace;
}

Face* HEdge::setFace(Face* f) {
	mFace = f;
	return mFace;
}

bool HEdge::flag() const {
	return mFlag;
}

bool HEdge::setFlag(bool b) {
	mFlag = b;
	return mFlag;
}

bool HEdge::isBoundary() const {
	return mBoundary;
}

bool HEdge::isValid() const {
	return mValid;
}

bool HEdge::setValid(bool b) {
	mValid = b;
	return mValid;
}
Vertex* HEdge::edgeVertex() const {
	Eigen::Vector3f new_vec(0,0,0);
	new_vec += mStart->position();
	new_vec += mNext->start()->position();
	new_vec += mFace->halfEdge()->start()->position();
	new_vec += mFace->halfEdge()->end()->position();
	new_vec /= 4.0;
	Vertex* v = new Vertex(new_vec);
	return v;
}
OneRingHEdge::OneRingHEdge(const Vertex* v) {
	if (v == nullptr) {
		mStart = nullptr;
		mNext = nullptr;
	} else {
		mStart = v->halfEdge();
		mNext = v->halfEdge();
	}
}

HEdge* OneRingHEdge::nextHEdge() {
	HEdge* ret = mNext;
	if (mNext != nullptr && mNext->prev()->twin() != mStart) {
		mNext = mNext->prev()->twin();
	} else {
		mNext = nullptr;
	}
	return ret;
}

OneRingVertex::OneRingVertex(const Vertex* v): ring(v) {
}

Vertex* OneRingVertex::nextVertex() {
	HEdge* he = ring.nextHEdge();
	return he != nullptr ? he->end() : nullptr;
}

Vertex::Vertex() : mHEdge(nullptr), mFlag(0) {
	mPosition = Eigen::Vector3f::Zero();
	mColor = VCOLOR_BLUE;
	mNormal = Eigen::Vector3f::Zero();
}

Vertex::Vertex(const Eigen::Vector3f& v): mPosition(v), mHEdge(nullptr), mFlag(0) {
	mColor = VCOLOR_BLUE;
	mNormal = Eigen::Vector3f::Zero();
}

Vertex::Vertex(float x, float y, float z): mHEdge(nullptr), mFlag(0) {
	mPosition = Eigen::Vector3f(x, y, z);
	mColor = VCOLOR_BLUE;
	mNormal = Eigen::Vector3f::Zero();
}


const Eigen::Vector3f& Vertex::position() const {
	return mPosition;
}

const Eigen::Vector3f& Vertex::setPosition(const Eigen::Vector3f& p) {
	mPosition = p;
	return mPosition;
}

const Eigen::Vector3f& Vertex::normal() const {
	return mNormal;
}

const Eigen::Vector3f& Vertex::setNormal(const Eigen::Vector3f& n) {
	mNormal = n;
	return mNormal;
}

const Eigen::Vector3f& Vertex::color() const {
	return mColor;
}

const Eigen::Vector3f& Vertex::setColor(const Eigen::Vector3f& c) {
	mColor = c;
	return mColor;
}

HEdge* Vertex::halfEdge() const {
	return mHEdge;
}

HEdge* Vertex::setHalfEdge(HEdge* he) {
	mHEdge = he;
	return mHEdge;
}

int Vertex::index() const {
	return mIndex;
}

int Vertex::setIndex(int i) {
	mIndex = i;
	return mIndex;
}

int Vertex::flag() const {
	return mFlag;
}

int Vertex::setFlag(int f) {
	mFlag = f;
	return mFlag;
}

bool Vertex::isValid() const {
	return mValid;
}

bool Vertex::setValid(bool b) {
	mValid = b;
	return mValid;
}

bool Vertex::isBoundary() const {
	OneRingHEdge ring(this);
	HEdge* curr = nullptr;
	while (curr = ring.nextHEdge()) {
		if (curr->isBoundary()) {
			return true;
		}
	}
	return false;
}

int Vertex::valence() const {
	int count = 0;
	OneRingVertex ring(this);
	Vertex* curr = nullptr;
	while (curr = ring.nextVertex()) {
		++count;
	}
	return count;
}

Face::Face() : mHEdge(nullptr), mValid(true) {
}

HEdge* Face::halfEdge() const {
	return mHEdge;
}
Vertex* Face::faceVertex() const {
	int num = 0;
	Eigen::Vector3f new_vec(0,0,0);
	HEdge* p = mHEdge;
	do
	{
		num++;
		new_vec += p->end()->position();
		p = p->next();
	}while( p!=mHEdge);
	new_vec /= num;
	Vertex* v = new Vertex(new_vec);
	return v;
}

HEdge* Face::setHalfEdge(HEdge* he) {
	mHEdge = he;
	return mHEdge;
}

bool Face::isBoundary() const {
	HEdge* curr = mHEdge;
	do {
		if (curr->twin()->isBoundary()) {
			return true;
		}
		curr = curr->next();
	} while (curr != mHEdge);
	return false;
}

bool Face::isValid() const {
	return mValid;
}

bool Face::setValid(bool b) {
	mValid = b;
	return mValid;
}

Mesh::Mesh() {
	mVertexPosFlag = true;
	mVertexNormalFlag = true;
	mVertexColorFlag = true;
}

Mesh::~Mesh() {
	clear();
}

const std::vector< HEdge* >& Mesh::edges() const {
	return mHEdgeList;
}

const std::vector< HEdge* >& Mesh::boundaryEdges() const {
	return mBHEdgeList;
}

const std::vector< Vertex* >& Mesh::vertices() const {
	return mVertexList;
}

const std::vector< Face* >& Mesh::faces() const {
	return mFaceList;
}


bool Mesh::isVertexPosDirty() const {
	return mVertexPosFlag;
}

void Mesh::setVertexPosDirty(bool b) {
	mVertexPosFlag = b;
}

bool Mesh::isVertexNormalDirty() const {
	return mVertexNormalFlag;
}

void Mesh::setVertexNormalDirty(bool b) {
	mVertexNormalFlag = b;
}

bool Mesh::isVertexColorDirty() const {
	return mVertexColorFlag;
}

void Mesh::setVertexColorDirty(bool b) {
	mVertexColorFlag = b;
}

bool Mesh::loadMeshFile(const std::string filename) {
	// Use libigl to parse the mesh file
	bool iglFlag = igl::read_triangle_mesh(filename, mVertexMat, mFaceMat);
	if (iglFlag) {
		clear();

		// Construct the half-edge data structure.
		int numVertices = mVertexMat.rows();
		int numFaces = mFaceMat.rows();

		// Fill in the vertex list
		for (int vidx = 0; vidx < numVertices; ++vidx) {
			mVertexList.push_back(new Vertex(mVertexMat(vidx, 0),
			                                 mVertexMat(vidx, 1),
			                                 mVertexMat(vidx, 2)));
		}
		// Fill in the face list
		for (int fidx = 0; fidx < numFaces; ++fidx) {
			addFace(mFaceMat(fidx, 0), mFaceMat(fidx, 1), mFaceMat(fidx, 2));
		}

		std::vector< HEdge* > hedgeList;
		for (int i = 0; i < mBHEdgeList.size(); ++i) {
			if (mBHEdgeList[i]->start()) {
				hedgeList.push_back(mBHEdgeList[i]);
			}
			// TODO
		}
		mBHEdgeList = hedgeList;

		for (int i = 0; i < mVertexList.size(); ++i) {
			mVertexList[i]->adjHEdges.clear();
			mVertexList[i]->setIndex(i);
			mVertexList[i]->setFlag(0);
		}
	} else {
		std::cout << __FUNCTION__ << ": mesh file loading failed!\n";
	}
	return iglFlag;
}

static void _setPrevNext(HEdge* e1, HEdge* e2) {
	e1->setNext(e2);
	e2->setPrev(e1);
}

static void _setTwin(HEdge* e1, HEdge* e2) {
	e1->setTwin(e2);
	e2->setTwin(e1);
}

static void _setFace(Face* f, HEdge* e) {
	f->setHalfEdge(e);
	e->setFace(f);
}

void Mesh::addFace(int v1, int v2, int v3) {
	Face* face = new Face();

	HEdge* hedge[3];
	HEdge* bhedge[3]; // Boundary half-edges
	Vertex* vert[3];

	for (int i = 0; i < 3; ++i) {
		hedge[i] = new HEdge();
		bhedge[i] = new HEdge(true);
	}
	vert[0] = mVertexList[v1];
	vert[1] = mVertexList[v2];
	vert[2] = mVertexList[v3];

	// Connect prev-next pointers
	for (int i = 0; i < 3; ++i) {
		_setPrevNext(hedge[i], hedge[(i + 1) % 3]);
		_setPrevNext(bhedge[i], bhedge[(i + 1) % 3]);
	}

	// Connect twin pointers
	_setTwin(hedge[0], bhedge[0]);
	_setTwin(hedge[1], bhedge[2]);
	_setTwin(hedge[2], bhedge[1]);

	// Connect start pointers for bhedge
	bhedge[0]->setStart(vert[1]);
	bhedge[1]->setStart(vert[0]);
	bhedge[2]->setStart(vert[2]);
	for (int i = 0; i < 3; ++i) {
		hedge[i]->setStart(vert[i]);
	}

	// Connect start pointers
	// Connect face-hedge pointers
	for (int i = 0; i < 3; ++i) {
		vert[i]->setHalfEdge(hedge[i]);
		vert[i]->adjHEdges.push_back(hedge[i]);
		_setFace(face, hedge[i]);
	}
	vert[0]->adjHEdges.push_back(bhedge[1]);
	vert[1]->adjHEdges.push_back(bhedge[0]);
	vert[2]->adjHEdges.push_back(bhedge[2]);

	// Merge boundary if needed
	for (int i = 0; i < 3; ++i) {
		Vertex* start = bhedge[i]->start();
		Vertex* end = bhedge[i]->end();

		for (int j = 0; j < end->adjHEdges.size(); ++j) {
			HEdge* curr = end->adjHEdges[j];
			if (curr->isBoundary() && curr->end() == start) {
				_setPrevNext(bhedge[i]->prev(), curr->next());
				_setPrevNext(curr->prev(), bhedge[i]->next());
				_setTwin(bhedge[i]->twin(), curr->twin());
				bhedge[i]->setStart(nullptr); // Mark as unused
				curr->setStart(nullptr); // Mark as unused
				break;
			}
		}
	}

	// Finally add hedges and faces to list
	for (int i = 0; i < 3; ++i) {
		mHEdgeList.push_back(hedge[i]);
		mBHEdgeList.push_back(bhedge[i]);
	}
	mFaceList.push_back(face);
}

Eigen::Vector3f Mesh::initBboxMin() const {
	return (mVertexMat.colwise().minCoeff()).transpose();
}

Eigen::Vector3f Mesh::initBboxMax() const {
	return (mVertexMat.colwise().maxCoeff()).transpose();
}

void Mesh::groupingVertexFlags() {
	// Init to 255
	for (Vertex* vert : mVertexList) {
		if (vert->flag() != 0) {
			vert->setFlag(255);
		}
	}
	// Group handles
	int id = 0;
	std::vector< Vertex* > tmpList;
	for (Vertex* vert : mVertexList) {
		if (vert->flag() == 255) {
			++id;
			vert->setFlag(id);

			// Do search
			tmpList.push_back(vert);
			while (!tmpList.empty()) {
				Vertex* v = tmpList.back();
				tmpList.pop_back();

				OneRingVertex orv = OneRingVertex(v);
				while (Vertex* v2 = orv.nextVertex()) {
					if (v2->flag() == 255) {
						v2->setFlag(id);
						tmpList.push_back(v2);
					}
				}
			}
		}
	}
}

void Mesh::clear() {
	for (int i = 0; i < mHEdgeList.size(); ++i) {
		delete mHEdgeList[i];
	}
	for (int i = 0; i < mBHEdgeList.size(); ++i) {
		delete mBHEdgeList[i];
	}
	for (int i = 0; i < mVertexList.size(); ++i) {
		delete mVertexList[i];
	}
	for (int i = 0; i < mFaceList.size(); ++i) {
		delete mFaceList[i];
	}

	mHEdgeList.clear();
	mBHEdgeList.clear();
	mVertexList.clear();
	mFaceList.clear();
}

std::vector< int > Mesh::collectMeshStats() {
	int V = 0; // # of vertices
	int E = 0; // # of half-edges
	int F = 0; // # of faces
	int B = 0; // # of boundary loops
	int C = 0; // # of connected components
	int G = 0; // # of genus

	/*====== Programming Assignment 0 ======*/

	/**********************************************/
	/*          Insert your code here.            */
	/**********************************************/
	/*
	/* Collect mesh information as listed above.
	/**********************************************/

	/*====== Programming Assignment 0 ======*/

	std::vector< int > stats;
	stats.push_back(V);
	stats.push_back(E);
	stats.push_back(F);
	stats.push_back(B);
	stats.push_back(C);
	stats.push_back(G);
	return stats;
}

int Mesh::countBoundaryLoops() {
	int count = 0;

	/*====== Programming Assignment 0 ======*/

	/**********************************************/
	/*          Insert your code here.            */
	/**********************************************/
	/*
	/* Helper function for Mesh::collectMeshStats()
	/**********************************************/

	/*====== Programming Assignment 0 ======*/

	return count;
}

int Mesh::countConnectedComponents() {
	int count = 0;

	/*====== Programming Assignment 0 ======*/

	/**********************************************/
	/*          Insert your code here.            */
	/**********************************************/
	/*
	/* Helper function for Mesh::collectMeshStats()
	/* Count the number of connected components of
	/* the mesh. (Hint: use a stack)
	/**********************************************/


	/*====== Programming Assignment 0 ======*/

	return count;
}

void Mesh::computeVertexNormals() {
	/*====== Programming Assignment 0 ======*/

	/**********************************************/
	/*          Insert your code here.            */
	/**********************************************/
	/*
	/* Compute per-vertex normal using neighboring
	/* facet information. (Hint: remember using a 
	/* weighting scheme. Plus, do you notice any
	/* disadvantages of your weighting scheme?)
	/**********************************************/

	/*====== Programming Assignment 0 ======*/
	
	// Notify mesh shaders
	setVertexNormalDirty(true);
}


void Mesh::umbrellaSmooth(bool cotangentWeights) {
	/*====== Programming Assignment 1 ======*/
	double lambda = 0.5;

	if (cotangentWeights) {
		/**********************************************/
		/*          Insert your code here.            */
		/**********************************************/
		/*
		/* Step 1: Implement the cotangent weighting 
		/* scheme for explicit mesh smoothing. 
		/*
		/* Hint:
		/* It is advised to double type to store the 
		/* weights to avoid numerical issues.
		/**********************************************/
		for(auto iter = mVertexList.begin(); iter != mVertexList.end(); iter++){
			const Eigen::Vector3f& pos = (*iter)->position();
			std::vector<HEdge*> e; // Edge list for the vertex
			// Add edges to the list.
			e.push_back((*iter)->halfEdge());
			while (e.back()->twin()->next()!=e[0]) {
				e.push_back(e.back()->twin()->next());
			}
			// Sum up the position of each neighbor according to the edge list.
			Eigen::Vector3f sumCotWVec(0, 0, 0);
			double cotw = 0.0;
			for(auto niter = e.begin(); niter != e.end(); niter++){
				double cotangent = (((*niter)->end()->position() - (*niter)->next()->end()->position()).dot( \
									(*niter)->next()->next()->end()->position()-(*niter)->next()->end()->position()) / \
									((*niter)->end()->position() - (*niter)->next()->end()->position()).cross( \
									(*niter)->next()->next()->end()->position()-(*niter)->next()->end()->position()).norm() + \
							((*niter)->twin()->end()->position() - (*niter)->twin()->next()->end()->position()).dot( \
									(*niter)->twin()->next()->next()->end()->position()-(*niter)->twin()->next()->end()->position()) / \
									((*niter)->twin()->end()->position() - (*niter)->twin()->next()->end()->position()).cross( \
									(*niter)->twin()->next()->next()->end()->position()-(*niter)->twin()->next()->end()->position()).norm()) / 2.0;
				cotw += cotangent;
				sumCotWVec += cotangent*(*niter)->end()->position();
			}
			// Calculate the laplician smoothing and update the new postion.
			Eigen::Vector3f laplacian = sumCotWVec / cotw - pos;
			(*iter)->setPosition(pos + lambda * laplacian);
		}
	} else {
		/**********************************************/
		/*          Insert your code here.            */
		/**********************************************/
		/*
		/* Step 2: Implement the uniform weighting 
		/* scheme for explicit mesh smoothing.
		/**********************************************/
		for(auto iter = mVertexList.begin(); iter != mVertexList.end(); iter++){
			const Eigen::Vector3f& pos = (*iter)->position();
			std::vector<HEdge*> e; // Edge list for the vertex
			// Add edges to the list.
			e.push_back((*iter)->halfEdge());
			while (e.back()->twin()->next()!=e[0]) {
				e.push_back(e.back()->twin()->next());
			}
			// Sum up the position of each neighbor according to the edge list.
			Eigen::Vector3f sumVec(0, 0, 0);
			for(auto niter = e.begin(); niter != e.end(); niter++) sumVec += (*niter)->end()->position();
			// Calculate the laplician smoothing and update the new postion.
			Eigen::Vector3f laplacian = sumVec / e.size() - pos;
			(*iter)->setPosition(pos + lambda * laplacian);
		}
	}

	/*====== Programming Assignment 1 ======*/

	computeVertexNormals();
	// Notify mesh shaders
	setVertexPosDirty(true);
}

void Mesh::implicitUmbrellaSmooth(bool cotangentWeights) {
	/*====== Programming Assignment 1 ======*/
	double lambda = 1;
	/* A sparse linear system Ax=b solver using the conjugate gradient method. */
	auto fnConjugateGradient = [](const Eigen::SparseMatrix< float >& A,
	                              const Eigen::VectorXf& b,
	                              int maxIterations,
	                              float errorTolerance,
	                              Eigen::VectorXf& x)
	{
		/**********************************************/
		/*          Insert your code here.            */
		/**********************************************/
		/*
		/* Params:
		/*  A: 
		/*  b: 
		/*  maxIterations:	Max number of iterations
		/*  errorTolerance: Error tolerance for the early stopping condition
		/*  x:				Stores the final solution, but should be initialized. 
		/**********************************************/
		/*
		/* Step 1: Implement the biconjugate gradient
		/* method.
		/* Hint: https://en.wikipedia.org/wiki/Biconjugate_gradient_method
		/**********************************************/
		// Follow the hint and implement the unpreconditioned version of the algorithm
		Eigen::VectorXf r = b - A * x;
		Eigen::VectorXf p = r;
        for(int k=0; k<maxIterations; ++k){
        	double alpha = r.norm() * r.norm() / (p.adjoint() * A * p);
        	Eigen::VectorXf xk = x + alpha * p;
        	if ((xk.norm()-x.norm()) < errorTolerance) break;
        	x = xk;
        	Eigen::VectorXf rk = r - alpha * A * p;
        	double beta = rk.norm() * rk.norm() / (r.norm() * r.norm());
        	r = rk;
        	p = r + beta * p;
        }
	};

	/* IMPORTANT:
	/* Please refer to the following link about the sparse matrix construction in Eigen. */
	/* http://eigen.tuxfamily.org/dox/group__TutorialSparse.html#title3 */

	if (cotangentWeights) {
		/**********************************************/
		/*          Insert your code here.            */
		/**********************************************/
		/*
		/* Step 2: Implement the cotangent weighting 
		/* scheme for implicit mesh smoothing. Use
		/* the above fnConjugateGradient for solving
		/* sparse linear systems.
		/*
		/* Hint:
		/* It is advised to double type to store the
		/* weights to avoid numerical issues.
		/**********************************************/
		Eigen::VectorXf b(3 * mVertexList.size());
		for (auto iter=mVertexList.begin(); iter!=mVertexList.end(); iter++) {
			for (int i=0; i<3; i++) {
				b[(*iter)->index() * 3 + i] = (*iter)->position()[i];
			}
		}
		typedef Eigen::Triplet<double> T;
		std::vector< T > tripletList;
		for(int i=0; i< 3 * mVertexList.size(); i++) tripletList.push_back(T(i,i,1 + lambda));
		for(auto iter = mVertexList.begin(); iter != mVertexList.end(); iter++){
			const Eigen::Vector3f& pos = (*iter)->position();
			std::vector<HEdge*> e; // Edge list for the vertex
			// Add edges to the list.
			e.push_back((*iter)->halfEdge());
			while (e.back()->twin()->next()!=e[0]) {
				e.push_back(e.back()->twin()->next());
			}
			std::vector<double> cotw;
			double sumcotw = 0.0;
			for(auto niter = e.begin(); niter != e.end(); niter++){
				double cotangent = (((*niter)->end()->position() - (*niter)->next()->end()->position()).dot( \
									(*niter)->next()->next()->end()->position()-(*niter)->next()->end()->position()) / \
									((*niter)->end()->position() - (*niter)->next()->end()->position()).cross( \
									(*niter)->next()->next()->end()->position()-(*niter)->next()->end()->position()).norm() + \
							((*niter)->twin()->end()->position() - (*niter)->twin()->next()->end()->position()).dot( \
									(*niter)->twin()->next()->next()->end()->position()-(*niter)->twin()->next()->end()->position()) / \
									((*niter)->twin()->end()->position() - (*niter)->twin()->next()->end()->position()).cross( \
									(*niter)->twin()->next()->next()->end()->position()-(*niter)->twin()->next()->end()->position()).norm()) / 2.0;
				cotw.push_back(cotangent);
				sumcotw += cotangent;
			}
			auto witer = cotw.begin();
			auto niter = e.begin();
			for(; niter != e.end(); niter++, witer++){
				for (int i=0; i<3; i++) {
					tripletList.push_back(T((*iter)->index() * 3 + i, (*niter)->end()->index() * 3 + i, -lambda * (*witer) / sumcotw));
				}
			}
		}
		//Solve the sparse linear problem.
		Eigen::SparseMatrix<float> A(3 * mVertexList.size(),3 * mVertexList.size());
		A.setFromTriplets(tripletList.begin(), tripletList.end());
		Eigen::VectorXf x(3 * mVertexList.size());
		for (int i=0; i<3*mVertexList.size(); i++) x[i] = 0.0;
		fnConjugateGradient(A, b, 15, 1e-7, x);
		// Calculate the laplician smoothing and update the new postion.
		for (auto iter=mVertexList.begin(); iter!=mVertexList.end(); iter++) {
			Eigen::Vector3f newpos;
			for (int i=0; i<3; i++) {
				newpos[i] = x[(*iter)->index() * 3 + i];
			}
			(*iter)->setPosition(newpos);
		}

	} else {
		/**********************************************/
		/*          Insert your code here.            */
		/**********************************************/
		/*
		/* Step 3: Implement the uniform weighting 
		/* scheme for implicit mesh smoothing. Use
		/* the above fnConjugateGradient for solving
		/* sparse linear systems.
		/**********************************************/
		Eigen::VectorXf b(3 * mVertexList.size());
		for (auto iter=mVertexList.begin(); iter!=mVertexList.end(); iter++) {
			for (int i=0; i<3; i++) {
				b[(*iter)->index() * 3 + i] = (*iter)->position()[i];
			}
		}
		typedef Eigen::Triplet<double> T;
		std::vector< T > tripletList;
		for(int i=0; i< 3 * mVertexList.size(); i++) tripletList.push_back(T(i,i,1 + lambda));
		for(auto iter = mVertexList.begin(); iter != mVertexList.end(); iter++){
			const Eigen::Vector3f& pos = (*iter)->position();
			std::vector<HEdge*> e; // Edge list for the vertex
			// Add edges to the list.
			e.push_back((*iter)->halfEdge());
			while (e.back()->twin()->next()!=e[0]) {
				e.push_back(e.back()->twin()->next());
			}
			for(auto niter = e.begin(); niter != e.end(); niter++){
				for (int i=0; i<3; i++) {
					tripletList.push_back(T((*iter)->index() * 3 + i, (*niter)->end()->index() * 3 + i, -lambda / e.size()));
				}
			}
		}
		//Solve the sparse linear problem.
		Eigen::SparseMatrix<float> A(3 * mVertexList.size(),3 * mVertexList.size());
		A.setFromTriplets(tripletList.begin(), tripletList.end());
		Eigen::VectorXf x(3 * mVertexList.size());
		for (int i=0; i<3*mVertexList.size(); i++) x[i] = 0.0;
		fnConjugateGradient(A, b, 15, 1e-7, x);
		// Calculate the laplician smoothing and update the new postion.
		for (auto iter=mVertexList.begin(); iter!=mVertexList.end(); iter++) {
			Eigen::Vector3f newpos;
			for (int i=0; i<3; i++) {
				newpos[i] = x[(*iter)->index() * 3 + i];
			}
			(*iter)->setPosition(newpos);
		}
	}

	/*====== Programming Assignment 1 ======*/

	computeVertexNormals();
	// Notify mesh shaders
	setVertexPosDirty(true);
}
