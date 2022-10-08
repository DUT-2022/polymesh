#include <iostream>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS.hxx>

#include "CrossFace.h"
#include "CrossFieldParams.h"
#include "create_msh_from_API.h"

#include "OCCD_Basic.h"

#include "PolyMesh.h"
#include "CrossfieldCommon.h"

using namespace std;
//using namespace Eigen;


void main(int argc, char** argv)
{
	string rootPath = "C:/Users/DUTuser/Desktop/";
	string geoName = "geo-curve";
	auto params = CrossFieldParams();
	params.bgMeshSize = 4;
	params.nIter = 1;
	params.partialBlockBool = true;
	cout << "Initialization..." << endl;

	// 1-1. Create background triangle mesh
	create_msh_from_API(rootPath, geoName);

	// 1-2. Read geometry file
	//暂时没有看明白，需要解决的问题较多
	TopologyUtils topo = read_step_file(rootPath + geoName + ".stp");

	//1-3 Read background tri mesh
	cout << 1 << endl;
	PolyMesh pm;
	pm.topo = topo;
	pm.loadMshFileTri(rootPath + geoName + ".msh");
	cout << 2 << endl;
	int i = 0;
	auto faces = topo.faces();
	for (TopoDS_Face face : faces)
	{
		
		list <MeshNode> nodesOnVert, nodesOnEdge, nodesInFace;
		list <MeshElement> faceElems;
		list <MeshEdge> faceMeshEdges;
		CrossFace cf(params);
		cf.getFaceNodes(pm, face, nodesOnVert, nodesOnEdge, nodesInFace);
		cf.getFaceElems(pm, i, faceElems);
		cf.getFaceMeshEdges(faceElems, faceMeshEdges);

		list <MeshNode> faceNodes;
		list <MeshNode> boundaryNodes;
		
		//可能有更好的方法来解决这个问题，[nodesOnVert,nodesOnEdge,nodesInFace]这是python程序，发挥你们的聪明才智吧，老学长太拉了，只能用笨方法。
		list <list<MeshNode>> faceNodesList;
		faceNodesList.push_back(nodesOnVert);
		faceNodesList.push_back(nodesOnEdge);
		faceNodesList.push_back(nodesInFace);
		
		faceNodes = flatten_list(faceNodesList);

		list <list<MeshNode>> boundaryNodesList;
		boundaryNodesList.push_back(nodesOnVert);
		boundaryNodesList.push_back(nodesOnEdge);

		boundaryNodes = flatten_list(boundaryNodesList);
		
		//cf.setNodeNormal(nodesOnVert, nodesOnEdge, nodesInFace, face);
	}
	
}
