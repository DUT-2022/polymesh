#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <Eigen/Dense>
#include <iterator>
#include <vector>
#include <list>
#include <time.h>

#include <conio.h>
#include "PolyMesh.h"
#include "Primitive3D.h"
#include "viewer.h"


#include <GL/glut.h>
#include <cstdio>
#include <cmath>


using namespace std;
using namespace Eigen;


MeshNode::MeshNode()
{

}


MeshElement::MeshElement()
{

}
MeshEdge::MeshEdge()
{

}

inline bool operator == (const MeshElement& p1, const MeshElement& p2)
{
	if (p2.mesh == p1.mesh && p1.ID == p2.ID && p1.meshEdges == p2.meshEdges && p1.startV == p2.startV && p1.normal == p2.normal)
		return true;
	return false;
}

inline bool operator == (const MeshEdge& p1, const MeshEdge& p2)
{
	if (p2.mesh == p1.mesh && p1.ID == p2.ID && p1.v1 == p2.v1 && p1.v2 == p2.v2 && p1.f1 == p2.f1 && p1.f2 == p2.f2)
		return true;
	return false;
}

inline bool operator == (const MeshNode& p1, const MeshNode& p2)
{
	if (p2.mesh == p1.mesh && p1.ID == p2.ID && p1.meshEdges == p2.meshEdges && p1.edgeIndex == p2.edgeIndex && p1.pos == p2.pos)
		return true;
	return false;
}

MeshNode::MeshNode(PolyMesh& M, Vector3d& P, int I)
{
	this->mesh = &M;
	this->pos = P;
	this->ID = I;
}

auto MeshNode::getNodeNeighbors()
{
	list<MeshNode*> ret;
	for (auto meshEdge : this->meshEdges)
	{
		ret.emplace_back(meshEdge->nodeAcross(*this));
	}
	return ret;
}

auto MeshNode::getAttachedElement()
{
	list<MeshElement*> ret;
	for (auto meshEdge : this->meshEdges)
	{
		if (meshEdge->f1 != NULL)
			ret.emplace_back(meshEdge->f1);
		if (meshEdge->f2 != NULL)
			ret.emplace_back(meshEdge->f2);
	}
	/*注意！！！，下方的sort函数会将原有顺序打乱，以便于unique函数会删除掉重复的元素*/
	ret.sort();
	ret.unique();
	return ret;
}

auto MeshNode::getOneRingArea()
{
	list<MeshElement*> elements = this->getAttachedElement();
	float ret = 0.0;
	for (auto element : elements)
	{
		ret = ret + element->getArea();
	}
	return ret;
}

auto MeshNode::getNormal()
{
	list<MeshElement*> elements = this->getAttachedElement();
	float totalArea = 0.0;
	Vector3d normal(0, 0, 0);
	for (auto element : elements)
	{
		auto w = element->getArea();
		totalArea = totalArea + w;
		normal = normal + w * element->getNormal();
	}
	if (totalArea == 0)
		return normal;
	Vector3d Normal = (1.0 / totalArea) * normal;
	return Normal;
}

void MeshNode::setPov(Vector3d& newPov)
{
	this->pos = newPov;
}

MeshElement::MeshElement(PolyMesh& M, int I)
{
	this->mesh = &M;
	this->ID = I;
}

void MeshElement::flipNormal()
{
	this->meshEdges.reverse();
	this->normal.clear();
}

auto MeshElement::getNodes()
{
	list<MeshNode*> ret;
	MeshNode* v = this->startV;
	for (auto meshEdge : this->meshEdges)
	{
		ret.emplace_back(v);
		v = meshEdge->nodeAcross(*v);
	}
	return ret;
}

auto MeshElement::getNodesPos()
{
	list<MeshNode*> ret = this->getNodes();
	list<Vector3d> pos;
	for (auto v : ret)
	{
		pos.emplace_back(v->pos);
	}
	return pos;
}

Vector3d MeshElement::getNormal()
{
	return getFaceNormal(this->getNodesPos());
}

float MeshElement::getArea()
{
	return getPolygonArea(this->getNodesPos());
}

auto MeshElement::getCentroid()
{
	list<Vector3d> poss = this->getNodesPos();
	Vector3d v(0, 0, 0);
	if (poss.size() == 0)
		return Vector3d((0, 0, 0));
	for (auto pos : poss)
	{
		v = v + pos;
	}
	Vector3d vv = v/ poss.size();
	return vv;
}

auto MeshElement::getPlane()
{
	return Plant3D(this->startV->pos, this->getNormal());
}

MeshEdge::MeshEdge(PolyMesh& M, MeshNode* v11, MeshNode* v22, int I)
{
	this->mesh = &M;
	this->ID = I;
	this->v1 = v11;
	this->v2 = v22;

}

MeshNode* MeshEdge::nodeAcross(MeshNode& startV)
{
	MeshNode* v = NULL;
	if (&startV == this->v1)
	{
		v = this->v2;
		return v;
	}
	else if (&startV == this->v2)
	{
		v = this->v1;
		return v;
	}
	else
		return v;
		cout << "警告（nodeAcross）：节点不是边的成员" << endl;
}

void MeshEdge::addElement(MeshElement& element)
{
	if (this->f1 == NULL)
	{
		this->f1 = &element;
	}
	else if (this->f2 == NULL)
	{
		this->f2 = &element;
	}
	else
		cout << "无法再添加单元，已有两个单元" << endl;
}

void MeshEdge::removeElement(MeshElement& element)
{
	if (this->f1 == &element)
		this->f1 = NULL;
	else if (this->f2 == &element)
		this->f2 = NULL;
	else
		cout << "无法删除指向从未属于边缘的元素的边缘指针" << endl;
}

auto MeshEdge::elementAcross(MeshElement& element)
{
	
	MeshElement* f = NULL;
	if (this->f1 == &element)
	{
		f = this->f2;
		return f;
	}
	else if (this->f2 == &element)
	{
		f = this->f1;
		return f;
	}
	else
	{
		return f;
		cout << "边不与该单元相邻" << endl;
	}

}

auto MeshEdge::getCenter()
{
	Vector3d v(0, 0, 0);
	v(0) = (this->v1->pos(0) + this->v2->pos(0)) / 2;
	v(1) = (this->v1->pos(1) + this->v2->pos(1)) / 2;
	v(2) = (this->v1->pos(2) + this->v2->pos(2)) / 2;
	return v;
}

auto MeshEdge::numAttachedElement()
{
	int ret = 0;
	if (this->f1 != NULL)
		ret++;
	if (this->f2 != NULL)
		ret++;
	return ret;
}

auto PolyMesh::getMeshEdge(MeshNode* v1, MeshNode* v2)
{
	MeshEdge* edge = NULL;
	for (auto e1 : v1->meshEdges)
	{
		for (auto e2 : v2->meshEdges)
		{
			if (e1 == e2)
			{
				edge = e1;
				return edge;
			}
		}
	}
	return edge;
}

auto PolyMesh::getBoundaryMeshEdges()
{
	for (auto edge = this->meshEdges.begin(); edge !=this->meshEdges.end();++edge)
	{
		int numElem = edge->numAttachedElement();
		if (numElem == 1)
			this->boundaryMeshEdges.emplace_back(&(*edge));
	}
	return this->boundaryMeshEdges;
}

auto PolyMesh::getBoundaryNodes()
{
	if (this->boundaryMeshEdges.size() == 0)
		this->getBoundaryMeshEdges();
	for (auto edge : this->boundaryMeshEdges)
	{
		int flag1 = 0;
		int flag2 = 0;
		for (auto node : this->boundaryNodes)
		{
			if (edge->v1 == node)
				flag1 = 1;
			if (edge->v2 == node)
				flag2 = 1;
		}
		if (flag1 == 0)
			this->boundaryNodes.emplace_back(edge->v1);
		if (flag2 == 0)
			this->boundaryNodes.emplace_back(edge->v2);
	}
	return this->boundaryNodes;
}

auto PolyMesh::getInternalNodes()
{
	if (this->boundaryNodes.size() == 0)
		this->getBoundaryNodes();
	for (auto node = this->nodes.begin(); node != this->nodes.end();++node)
	{
		int flag = 0;
		for (auto node1 : this->boundaryNodes)
		{
			if (&(* node) == node1)
				flag = 1;
		}
		for (auto node2 : this->internalNodes)
		{
			if (&(*node) == node2)
				flag = 1;
		}
		if (flag == 0)
			this->internalNodes.emplace_back(&(*node));
	}
	return this->internalNodes;
}

void PolyMesh::getOrderedBoundaryMeshEdges()
{
	if (this->boundaryMeshEdges.size() == 0)
		this->getBoundaryMeshEdges();
}


void PolyMesh::addNode(Vector3d& Pos)
{	
	MeshNode node(*this, Pos, this->nodes.size());
	this->nodes.emplace_back(node);
}

MeshEdge* PolyMesh::addMeshEdge(MeshNode* v1, MeshNode* v2)
{
	
	MeshEdge edge(*this,v1, v2, this->edgeIndex);
	this->meshEdges.emplace_back(edge);
	MeshEdge* edge1 = &this->meshEdges.back();
	v1->meshEdges.emplace_back(edge1);
	v2->meshEdges.emplace_back(edge1);
	this->edgeIndex++;
	
	return edge1;
}

void PolyMesh::addElement(vector<MeshNode*>& verts)
{
	MeshNode* v1;
	MeshNode* v2;
	MeshEdge* edge;
	MeshElement element(*this,this->elements.size());
	element.startV = verts[0];
	this->elements.emplace_back(element);
	for (int i = 0; i < verts.size(); i++)
	{
		v1 = verts[i];
		v2 = verts[(i + 1) % verts.size()];
		if (this->getMeshEdge(v1, v2) == NULL)
		{
			edge = this->addMeshEdge(v1, v2);
			this->elements.back().meshEdges.emplace_back(edge);
			edge->addElement(this->elements.back());
		}

		else
		{
			edge = this->getMeshEdge(v1, v2);
			this->elements.back().meshEdges.emplace_back(edge);
			edge->addElement(this->elements.back());
		}
			
	}

}

void PolyMesh::removeElement(MeshElement& element)
{
	this->elements.remove(element);
	for (auto edge : element.meshEdges)
	{
		edge->removeElement(element);
	}

}

void PolyMesh::removeMeshEdge(MeshEdge& edge)
{
	this->meshEdges.remove(edge);
	edge.v1->meshEdges.remove(&edge);
	edge.v2->meshEdges.remove(&edge);

}

void PolyMesh::removeNode(MeshNode& node)
{
	this->nodes.remove(node);
}

void PolyMesh::loadOffFile(string filename)
{
	int nNodes = 0;
	int nElements = 0;
	int nEdges = 0;
	int node = 0;
	int element = 0;
	ifstream infile;
	infile.open(filename, ios::in);
	//infile.open("C:/Users/26065/Desktop/ccccc.off", ios::in);
	if (!infile.is_open())
	{
		cout << "文件打开失败" << endl;
		exit(0);
	}
	string buf;
	while (getline(infile, buf))
	{
		istringstream in(buf);
		vector<string> v;
		string t;
		while (in >> t)
		{
			v.push_back(t);
		}
		int size = v.size();
		vector<double> d;
		for (int i = 0; i < v.size(); i++)
		{
			
			stringstream ss;
			double h;
			ss <<  v[i];
			ss  >>  h;
			ss.clear();
			d.push_back(h);
		}
		if (nNodes == 0)
		{
			if (v[0] == "OFF")
			{
				if (v.size() > 2)
				{
					nNodes = int(d[1]);
					nElements = int(d[2]);
					nEdges = int(d[3]);
					this->nNodes = nNodes;
					this->nElements = nElements;
				}
			}
			else
			{
				nNodes = int(d[0]);
				nElements = int(d[1]);
				nEdges = int(d[2]);
				this->nNodes = nNodes;
				this->nElements = nElements;
			}
		}

		else if (node < nNodes)
		{
			Vector3d P(d[0], d[1], d[2]);
			this->addNode(P);
			node++;	
				
		}
		else if (element < nElements && node >= nNodes)
		{	
			vector<MeshNode*> verts;
			for(int i = 0; i < d[0];i++)
			{
				auto iter = this->nodes.begin();
				advance(iter, d[i+1]);
				verts.push_back(&*iter);
			}
			this->addElement(verts);
			element++;
		}
		
	}
	infile.close();
	return ;
}

void PolyMesh::loadMshFileTri(string filename)
{
	int nNodes = 0;
	int nElements = 0;
	int nEdges = 0;
	int node = 0;
	int element = 0;
	int flag_node = 0;
	int flag_element = 0;
	ifstream infile;
	infile.open(filename, ios::in);
	if (!infile.is_open())
	{
		cout << "文件打开失败" << endl;
		exit(0);
	}
	else
		cout << "文件成功打开" << endl;
	string buf;
	while (getline(infile, buf))
	{
		istringstream in(buf);
		vector<string> v;
		string t;
		while (in >> t)
		{
			v.push_back(t);
		}
		int size = v.size();
		vector<double> d;
		for (int i = 0; i < v.size(); i++)
		{
			stringstream ss;
			double h;
			ss << setprecision(16) << v[i];
			ss >> setprecision(16) >> h;
			ss.clear();
			d.push_back(h);
		}
		if (flag_node == 0)
		{
			if (v[0] == "$ParametricNodes" || v[0] == "$Nodes")
				flag_node = 1;
		}
		else if (flag_node == 1)
		{
			flag_node = 2;
		}
		else if (flag_node == 2 && v.size() > 5)
		{
			int geoIndex = int(d[5]) - 1;
			Vector3d P(d[1], d[2], d[3]);
			this->addNode(P);

			this->nodes.back().geoIndex = geoIndex;
			
			if (v[4] == "0")
			{
				this->nodes.back().geoType = 0;
				this->vertIndex_nodeID.emplace_back(int(d[0]) - 1);//There may be a problem.
			}
			else if (v[4] == "1")
			{
				this->nodes.back().geoType = 1;
				this->nodes.back().geoParm[0] = float(d[6]);
				
				if (this->edgeIndex_nodeIDs.size() <= geoIndex)
				{
					
					this->edgeIndex_nodeIDs.emplace_back(list<int>(int(d[0]) - 1));//There may be a problem.
					//cout << this->faceIndex_nodeIDs.back().back();
				}
				else
				{
					
					auto iter = edgeIndex_nodeIDs.begin();
					advance(iter, geoIndex);
					iter->emplace_back(int(d[0]) - 1);
				}
			}
			else if (v[4] == "2")
			{
				this->nodes.back().geoType = 2;
				this->nodes.back().geoParm[0] = float(d[6]);
				this->nodes.back().geoParm[1] = float(d[7]);
				if (this->faceIndex_nodeIDs.size() <= geoIndex)
				{
					
					this->faceIndex_nodeIDs.emplace_back(list<int>(int(d[0]) - 1));//There may be a problem.
					
				}
				else
				{
					auto iter = faceIndex_nodeIDs.begin();
					advance(iter, geoIndex);
					iter->emplace_back(int(d[0]) - 1);
				}
			}
			else if (v[4] == "3")
			{
				cout << "Please check if this is a surface mesh!!" << endl;
				break;
			}
		}
		else if (v[0] == "$EndParametricNodes")
		{
			flag_node = 3;
		}
		if (flag_element == 0)
		{
			if (v[0] == "$Elements")
			flag_element = 1;
		}
		else if (flag_element == 1)
		{
			flag_element = 2;
		}
		else if (flag_element == 2 && v.size() > 5)
		{
			if (v[1] != "2")
				continue;	
			else
			{	
				vector<MeshNode*> verts;
				for (int i = 2; i > -1; i--)
				{
					auto iter = this->nodes.begin();
					advance(iter, int(d[d.size()- i - 1]) - 1);
					verts.push_back(&*iter);
				}
				this->addElement(verts);
				if (this->faceIndex_elementID.size() <= int(d[4]) - 1)
				{
					this->faceIndex_nodeIDs.emplace_back(list<int>(this->elements.size()));//There may be a problem.
				}
				else
				{
					auto iter = faceIndex_nodeIDs.begin();
					advance(iter, int(d[4]) - 1);
					iter->emplace_back(this->elements.size());
				}

			}
		}
		else if (v[0] == "$EndParametricNodes")
		{
			flag_element = 3;
			infile.close();
		}
	}

}

void PolyMesh::saveMshFileQuad(string filename)
{
	ofstream OutFile(filename + "test.msh");
	OutFile << "$MeshFormat\n";
	OutFile << "2.2 0 8\n"; 
	OutFile << "$EndMeshFormat\n";
	OutFile << "$Nodes\n";
	OutFile << to_string(this->nodes.size()) + "\n";
	for (auto node : this->nodes)
	{
		OutFile << to_string(node.ID + 1);
		for (auto p : node.pos)
		{
			string h;
			stringstream ss;
			ss << setprecision(16) << p;
			ss >> h;
			ss.clear();
			OutFile << " " + h;
		}
		OutFile<< "\n";
	}
		
	OutFile << "$EndNodes\n";
	OutFile << "$Elements\n";
	OutFile << to_string(this->elements.size()) + "\n";
	for (auto element : this->elements)
	{
		auto nodes = element.getNodes();
		OutFile << to_string(element.ID + 1) + " 3 2 0 2";
		for (auto node : nodes)
			OutFile << " " + to_string(node->ID + 1);
		OutFile << "\n";
	}
	OutFile << "$EndElements\n";
	OutFile.close();
}

void PolyMesh::saveMshFileTri(string filename)
{
	
	ofstream OutFile(filename + "test1.msh");
	OutFile << "$MeshFormat\n";
	OutFile << "2.2 0 8\n";
	OutFile << "$EndMeshFormat\n";
	OutFile << "$Nodes\n";

	OutFile << to_string(this->nodes.size()) + "\n";
	for (auto node : this->nodes)
		OutFile << to_string(node.ID + 1) + " " + to_string(node.pos[0]) + " " + to_string(node.pos[1]) + " " + to_string(node.pos[2]) + "\n";
	OutFile << "$EndNodes\n";

	OutFile << "$Elements\n";
	OutFile << to_string(this->elements.size()) + "\n";
	for (auto element : this->elements)
	{
		auto nodes = element.getNodes();
		OutFile << to_string(element.ID + 1) + " 2 2 0 1";
		for (auto node : nodes)
			OutFile << " " + to_string(node->ID + 1);
		OutFile << "\n";
	}
	OutFile << "$EndElements\n";
	OutFile.close();
	
}
/*

void main(int argc, char** argv)
{

	PolyMesh M;

	string filename = "C:/Users/26065/Desktop/Shared/1-code/decomposition_fem1_mesh1.off";
	//string filename = "C:/Users/26065/Desktop/test2.msh";
	string filename1 = "C:/Users/26065/Desktop/";

	M.loadOffFile(filename);
	//M.loadMshFileTri(filename);
	M.saveMshFileQuad(filename1);
	//M.saveMshFileTri(filename1);


	initOff (argc, argv, filename);



	

	
	return;
}
*/