#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS.hxx>

#include <TopExp_Explorer.hxx>
#include <BRep_Tool.hxx>
#include "TopologyUtils.h"
#include <tuple>
#include <vector>
using namespace std;


TopologyUtils::TopologyUtils()
{
	
}

TopologyUtils::TopologyUtils(TopoDS_Shape Topo)
{
	this->topo = Topo;
}

// Get faces
std::vector<TopoDS_Face> TopologyUtils::faces()
{
	std::vector<TopoDS_Face>  faces;
	faces.resize(0);
	TopExp_Explorer exp;
	for (exp.Init(this->topo, TopAbs_FACE); exp.More(); exp.Next())
	{
		faces.push_back(TopoDS::Face(exp.Current()));
	}

	return faces;
}



// Get edges
std::vector<TopoDS_Edge> TopologyUtils::edges()
{
	std::vector<TopoDS_Edge>  edges;
	edges.resize(0);
	TopExp_Explorer exp;
	for (exp.Init(this->topo, TopAbs_EDGE); exp.More(); exp.Next())
	{
		edges.push_back(TopoDS::Edge(exp.Current()));
	}

	return edges;
}
// Get vertices
std::vector<TopoDS_Vertex> TopologyUtils::vertices()
{
	std::vector<TopoDS_Vertex>  vertices;
	vertices.resize(0);
	TopExp_Explorer exp;
	for (exp.Init(this->topo, TopAbs_VERTEX); exp.More(); exp.Next())
	{
		vertices.push_back(TopoDS::Vertex(exp.Current()));
	}

	return vertices;
}