#include <iostream>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include "CrossFieldParams.h"
#include "CrossFace.h"
#include "PolyMesh.h"

using namespace std;

CrossFace::CrossFace(CrossFieldParams &params)
{
	this->crossFieldParams = params;
}

void CrossFace::getFaceElems(PolyMesh &pm, int faceIndex, list<MeshElement> &faceElems)
{
	//Get all the elements in a face

	if (pm.topo.faces().size() == 1)//Only one face
		faceElems = pm.elements;

}

void CrossFace::getFaceMeshEdges(list <MeshElement> &faceElems, list <MeshEdge> &faceMeshEdges)
{
	//Get all the mesh edges of a face
	for (auto elem : faceElems)
		for (auto edge : elem.meshEdges)
			faceMeshEdges.push_back(*edge);

}

void CrossFace::getFaceNodes(PolyMesh &pm, TopoDS_Face &face, list <MeshNode> &nodesOnVert, list <MeshNode> &nodesOnEdge, list <MeshNode> &nodesInFace)
{
	//Get nodes lie on the input face.
	if (pm.topo.faces().size() == 1)
	{
		for (auto node : pm.nodes)
		{
			if (node.geoType == 0)
				nodesOnVert.push_back(node);
			else if (node.geoType == 1)
				nodesOnEdge.push_back(node);
			else if (node.geoType == 2)
				nodesInFace.push_back(node);
		}
	}
}
/*
void setNodeNormal(list <MeshNode>nodesOnVert, list <MeshNode>nodesOnEdge, list <MeshNode>nodesInFace, TopoDS_Face face)
{
	//Calculate face normals at nodes. Note for nodes lies on vertex or
	//edges, its face normal value will be overrided when doing the calculation
	//for the next face.

	int typea = OCCB.ask_face_type(face);
	Vector3d normal;
	char uv;
	if (typea == 0)
	{
		if (nodesOnVert.size() != 0)
		{
			uv = OCCB.ask_point_uv2((nodesOnVert.front()).getPos(), face);
			normal = OCCB.ask_point_normal_face(uv, face);
		}
		else if (nodesOnEdge.size() != 0)
		{
			uv = OCCB.ask_point_uv2((nodesOnEdge.front()).getPos(), face);
			normal = OCCB.ask_point_normal_face(uv, face);
		}
		else if (nodesInFace.size() != 0)
		{
			Vector3d uv;
			if ((nodesInFace.front()).geoParm)
				uv = OCCB.ask_point_uv2((nodesInFace.front()).getPos(), face);

			else
				uv = (nodesInFace.front()).geoParm;
			normal = OCCB.ask_point_normal_face(uv, face);
		}
		else
			cout << "No nodes on the face, please check..." << endl;

		for (auto node : nodesOnVert)
			(node.faceNormal) = normal;
		for (auto node : nodesOnEdge)
			(node.faceNormal) = normal;
		for (auto node : nodesInFace)
			(node.faceNormal) = normal;
	}
	else

	{	//For nodes on vertex
		for (auto node : nodesOnVert)
		{
			uv = OCCB.ask_point_uv2(node.getPos(), face);
			normal = OCCB.ask_point_normal_face(uv, face);
			node.faceNormal = normal;
			//#print(normal);
		}

		//For nodes on edge
		for (auto node : nodesOnEdge)
		{
			uv = OCCB.ask_point_uv2(node.getPos(), face);
			normal = OCCB.ask_point_normal_face(uv, face);
			node.faceNormal = normal;
			//#print(normal);

		}

		//For nodes in face
		for (auto node : nodesInFace)
		{
			if (node.geoParm == [])
				uv = OCCB.ask_point_uv2(node.getPos(), face);
			else
				uv = node.geoParm;
			normal = OCCB.ask_point_normal_face(uv, face);
			node.faceNormal = normal;
		}
	}
	


}*/


                            