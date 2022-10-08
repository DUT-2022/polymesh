#pragma once
#include <tuple>
#include <vector>
#include <TopoDS_Shape.hxx>
#include <TopoDS.hxx>
using namespace std;
class TopologyUtils
{
    public:
        TopoDS_Shape topo;
		TopologyUtils();
        TopologyUtils(TopoDS_Shape);
        std::vector<TopoDS_Face> faces();
        // Get edges
        std::vector<TopoDS_Edge> edges();
        // Get vertices
        std::vector<TopoDS_Vertex> vertices();

};

