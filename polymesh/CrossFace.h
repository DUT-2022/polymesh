#pragma once
#include <iostream>
#include "CrossFieldParams.h"
#include "PolyMesh.h"
using namespace std;

class CrossFace
{
public:
	//Create a cross field on a face and use the cross field to identify singularities
	//and create split lines

	CrossFieldParams crossFieldParams;
	
	CrossFace(CrossFieldParams &params);
	void getFaceElems(PolyMesh &pm, int faceIndex, list<MeshElement> &faceElems);
	void getFaceMeshEdges(list <MeshElement> &faceElems, list <MeshEdge> &faceMeshEdges);
	void getFaceNodes(PolyMesh &pm, TopoDS_Face &face, list <MeshNode> &nodesOnVert, list <MeshNode> &nodesOnEdge, list <MeshNode> &nodesInFace);
	//void setNodeNormal(list <MeshNode>nodesOnVert, list <MeshNode>nodesOnEdge, list <MeshNode>nodesInFace, TopoDS_Face face)
};
