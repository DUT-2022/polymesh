#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <list>;
#include "TopologyUtils.h"

using namespace std;
using namespace Eigen;


class MeshNode;
class MeshElement;
class MeshEdge;
class PolyMesh
{
public:
	int nElements = 0;
	int nNodes = 0;
	int edgeIndex = 0;
	
	list<MeshNode*> boundaryNodes;
	list<MeshNode*> internalNodes;
	list<MeshEdge*> boundaryMeshEdges;
	list<MeshNode> nodes;
	list<MeshEdge> meshEdges;
	list<MeshElement> elements;
	list<int> vertIndex_nodeID;
	list<list<int>> edgeIndex_nodeIDs;
	list<list<int>> faceIndex_nodeIDs;
	list<list<int>> faceIndex_elementID;
	TopologyUtils topo;

	void loadOffFile(string);
	void loadMshFileTri(string);
	void saveMshFileQuad(string);
	void saveMshFileTri(string);
	auto getMeshEdge(MeshNode*, MeshNode*);
	void addNode(Vector3d&);
	MeshEdge* addMeshEdge(MeshNode* , MeshNode* );
	void addElement(vector<MeshNode*>&);
	void removeElement(MeshElement&);
	void removeMeshEdge(MeshEdge&);
	void removeNode(MeshNode&);
	auto getBoundaryMeshEdges();
	auto getBoundaryNodes();
	auto getInternalNodes();
	void getOrderedBoundaryMeshEdges();


};

class MeshNode
{
public:
	PolyMesh* mesh = NULL;
	int ID = 0;
	int edgeIndex = 0;
	Vector3d pos;
	list<MeshEdge*> meshEdges;
	int geoIndex = 0;
	int geoType = 0;
	float geoParm[2] = { 0.0,0.0 };

	auto getNodeNeighbors();
	auto getAttachedElement();
	auto getOneRingArea();
	auto getNormal();
	void setPov(Vector3d&);
	MeshNode(PolyMesh&, Vector3d&, int);
	MeshNode();
};

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


class MeshEdge
{
public:
	PolyMesh* mesh = NULL;
	int ID = 0;
	MeshNode* v1 = NULL;
	MeshNode* v2 = NULL;
	MeshElement* f1 = NULL;
	MeshElement* f2 = NULL;

	MeshNode* nodeAcross(MeshNode&);
	void addElement(MeshElement&);
	void removeElement(MeshElement&);
	auto elementAcross(MeshElement&);
	auto getCenter();
	auto numAttachedElement();

	MeshEdge(PolyMesh&, MeshNode* v1, MeshNode* v2, int I);
	MeshEdge();

};


