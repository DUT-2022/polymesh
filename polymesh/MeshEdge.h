#pragma once
#ifndef MeshEdge_H_
#define MeshEdge_H_
#include <Eigen/Dense>
#include <iostream>


using namespace std;
using namespace Eigen;


class MeshEdge
{
public:
	PolyMesh* mesh = NULL;
	int ID = 0;
	MeshNode* v1 = NULL;
	MeshNode* v2 = NULL;
	MeshElement* f1 = NULL;
	MeshElement* f2 = NULL;

	MeshNode* nodeAcross(MeshNode& );
	void addElement(MeshElement&);
	void removeElement(MeshElement&);
	auto elementAcross(MeshElement&);
	auto getCenter();
	auto numAttachedElement();

	MeshEdge(PolyMesh&, MeshNode* v1, MeshNode* v2, int I);
	MeshEdge();

};
MeshEdge::MeshEdge()
{

}

#endif