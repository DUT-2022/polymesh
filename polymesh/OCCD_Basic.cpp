#include <TopoDS_Shape.hxx>
#include <STEPControl_Reader.hxx>
#include "TopologyUtils.h"

#include <iostream>
using namespace std;

TopologyUtils read_step_file(const std::string& filename)
{
	STEPControl_Reader    reader;
	IFSelect_ReturnStatus stat;
	stat = reader.ReadFile(filename.c_str());
	//AssertThrow(stat == IFSelect_RetDone, ExcMessage("Error in reading file!"));

	auto    failsonly = true;
	IFSelect_PrintCount mode = IFSelect_ItemsByEntity;
	reader.PrintCheckLoad(failsonly, mode);
	Standard_Integer nRoots = reader.TransferRoots();
	// selects all IGES entities (including non visible ones) in the
	// file and puts them into a list called MyList,

	//AssertThrow(nRoots > 0, ExcMessage("Read nothing from file."));


	TopoDS_Shape sh = reader.OneShape();

	TopologyUtils topo(sh);
	
	return topo; // this is the actual translation
}