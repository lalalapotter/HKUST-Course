#include "deformer.h"
#include <iostream>

Deformer::Deformer() : mMesh(nullptr),
                       mCholeskySolver(nullptr) {
}

Deformer::~Deformer() {
	clear();
}

void Deformer::clear() {
	if (mCholeskySolver) {
		delete mCholeskySolver;
	}
	mCholeskySolver = nullptr;
	mRoiList.clear();
}

void Deformer::setMesh(Mesh* mesh) {
	mMesh = mesh;
	clear();
	// Record the handle vertices
	for (Vertex* vert : mMesh->vertices()) {
		if (vert->flag() > 0 || vert->isBoundary()) {
			mRoiList.push_back(vert);
		}
	}
	// Build system matrix for deformation
	buildSystemMat();
}


void Deformer::buildSystemMat() {
	/*====== Programming Assignment 2 ======*/

	/**********************************************/
	/*          Insert your code here.            */
	/**********************************************/
	/*
	/* Build the matrix of the linear system for 
	/* deformation and do factorization, in order
	/* to reuse and speed up in Deformer::deform().
	/* Handle vertices are maked by Vertex::flag() > 0
	/* Movements of the specified handle are already
	/* recorded in Vertex::position()
	/**********************************************/
	mVertexList = mMesh->vertices();

	Eigen::SparseMatrix< double > systemMat(mVertexList.size(), mVertexList.size());

	/*====== Programming Assignment 2 ======*/
	b = new Eigen::MatrixX3d (mVertexList.size()+mRoiList.size()+1, 3);
	Eigen::MatrixX3d b_original(mVertexList.size(), 3);
	for (auto iter=mVertexList.begin(); iter!=mVertexList.end(); iter++) {
		for (int i=0; i<3; i++) {
			b_original((*iter)->index(), i) = (*iter)->position()[i];
		}
	}
	typedef Eigen::Triplet<double> T;
	std::vector< T > tripletList;
	for(int i=0; i<mVertexList.size(); i++) tripletList.push_back(T(i,i,1));
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
			tripletList.push_back(T((*iter)->index(), (*niter)->end()->index(), -(*witer) / sumcotw));
		}
	}
	//Solve the sparse linear problem.
	A = new Eigen::SparseMatrix<double> (mVertexList.size()+mRoiList.size()+1,mVertexList.size());
	tripletList.push_back(T(mVertexList.size() , mVertexList[0]->index(), 1.0));
	int count = 0;
	for (auto iter=mRoiList.begin(); iter!=mRoiList.end(); iter++) {
		tripletList.push_back(T(mVertexList.size() + count + 1, (*iter)->index(), 1.0));
		count++;
	}	
	A->setFromTriplets(tripletList.begin(), tripletList.end());
	systemMat = A->transpose() * (*A);

	// Compute Laplacian coordinates
	*b = (*A) * b_original;
	// Please refer to the following link for the usage of sparse linear system solvers in Eigen
	// https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html

	// Do factorization
	if (systemMat.nonZeros() > 0) {
		mCholeskySolver = new Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > >();
		mCholeskySolver->compute(systemMat);
		if (mCholeskySolver->info() != Eigen::Success) {
			// Decomposition failed
			std::cout << "Sparse decomposition failed\n";
		} else {
			std::cout << "Sparse decomposition succeeded\n";
		}
	}
}

void Deformer::deform() {
	if (mCholeskySolver == nullptr) {
		return;
	}

	/*====== Programming Assignment 2 ======*/

	/**********************************************/
	/*          Insert your code here.            */
	/**********************************************/
	/*
	/* This is the place where the editing techniques 
	/* take place.
	/* Solve for the new vertex positions after the 
	/* specified handles move using the factorized
	/* matrix from Deformer::buildSystemMat(), i.e.,
	/* mCholeskySolver defined in deformer.h
	/**********************************************/
	// Please refer to the following link for the usage of sparse linear system solvers in Eigen
	// https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html

	// 
	mVertexList = mMesh->vertices();
	for (int i=0; i<3; i++) {
		(*b)(mVertexList.size(), i) = mVertexList[0]->position()[i];
	}
	int count = 0;
	for (auto it=mRoiList.begin(); it!=mRoiList.end(); it++) {
		Vertex* curVertex = *it;
		for (int i=0; i<3; i++) {
			(*b)(mVertexList.size() + count + 1, i) = curVertex->position()[i];
		}
		count++;
	}
	// Solve the matrix
	Eigen::MatrixX3d x(mVertexList.size(), 3);
	x = mCholeskySolver->solve(A->transpose() * (*b));
	// Set new position
	for (auto iter=mVertexList.begin(); iter!=mVertexList.end(); iter++) {
			Eigen::Vector3f newpos;
			for (int i=0; i<3; i++) {
				newpos[i] = x((*iter)->index(), i);
			}
			(*iter)->setPosition(newpos);
	}
	/*====== Programming Assignment 2 ======*/
}
