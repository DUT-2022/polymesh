#pragma once
#ifndef MeshNode_H_
#define MeshNode_H_
#include <Eigen/Dense>
#include <iostream>
#include <list>

using namespace std;
using namespace Eigen;




class MeshNode
{
public:
	PolyMesh* mesh = NULL;
	int ID = 0;
	int edgeIndex = 0;
	Vector3d pos ;
	list<MeshEdge*> meshEdges;
	int geoIndex = 0;
	int geoType = 0;
	float geoParm[2] = {0.0,0.0};

	auto getNodeNeighbors();
	auto getAttachedElement();
	auto getOneRingArea();
	auto getNormal();
	void setPov(Vector3d&);
	MeshNode( PolyMesh&, Vector3d& , int);
	MeshNode();
};
MeshNode::MeshNode()
{

}
#endif