#pragma once
#ifndef MeshElement_H_
#define MeshElement_H_
#include <Eigen/Dense>
#include <iostream>


using namespace std;
using namespace Eigen;



class MeshElement
{
public:
	PolyMesh* mesh = NULL;
	int ID = 0;
	list<MeshEdge*> meshEdges;
	MeshNode* startV = NULL;
	list<int> normal;
	

	MeshElement(PolyMesh&, int);
	MeshElement();
	void flipNormal();
	auto getNodes();
	auto getNodesPos();
	Vector3d getNormal();
	float getArea();
	auto getCentroid();
	auto getPlane();
};

MeshElement::MeshElement()
{

}
#endif