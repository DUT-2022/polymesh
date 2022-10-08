#include <iostream>
#include <string>
#include <stdio.h>
#include <io.h>
#include <cstring>
#include <fstream>
#include <math.h>
#include "create_msh_from_API.h"
#include <gmsh.h_cwrap>

using namespace std;

void create_msh_from_API(string path, string modelName)
{
	using namespace gmsh;

	int lc = 5;
	
	initialize();
	//open(path + modelName + ".stp");


	//gmsh.option.setNumber("General.Verbosity", 99)
	option::setNumber("Mesh.Algorithm", 1);
	//gmsh.option.setNumber('Mesh.MeshSizeMin', lc)
	option::setNumber("Mesh.MeshSizeMax", lc);
	//gmsh.option.setNumber("Mesh.MeshSizeMax", lc)
	option::setNumber("Mesh.MeshSizeFromCurvature", 100);
	option::setNumber("Mesh.SaveParametric", 1);
	//gmsh.option.setNumber("Mesh.Format", 1);
	option::setNumber("Mesh.MshFileVersion", 2.2);
	option::setNumber("Mesh.SaveAll", 1);
	option::setNumber("Mesh.MeshSizeExtendFromBoundary", 1);

	//Import BREP, STEP or IGES shapes from the file in the OpenCASCADE CAD representation
	std::vector<std::pair<int, int>> dimTags;
	model::occ::importShapes(path + modelName + ".stp", dimTags);
	std::cout << dimTags.size() << std::endl;
	for (std::size_t i = 0; i < dimTags.size(); i++)
	{
		std::cout << "(" << dimTags[i].first << ", " << dimTags[i].second << ")" << std::endl;

	}

	std::vector<std::pair<int, int> > entities;
	gmsh::model::getEntities(entities);

	for (auto e : entities) {
		// Dimension and tag of the entity:
		int dim = e.first, tag = e.second;
		std::cout << tag << std::endl;
	}
	model::occ::synchronize();

	//gmsh.model.geo.synchronize();
	model::mesh::generate(2);


	write(path + modelName + ".msh");

	//fltk::run();
	finalize();
}


