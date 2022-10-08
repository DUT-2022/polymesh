#include <iostream>
#include <list>
#include <Eigen/Dense>
#include "Primitive3D.h"

using namespace std;
using namespace Eigen;
double EPS = 0.000000000001;

auto normalizeVec(Vector3d& ret)
{
	ret.normalize();
}

Vector3d getFaceNormal(list<Vector3d> pos)
{
	auto pos0 = pos.begin();
	
	for (int i = 2; i < pos.size(); ++i)
	{
		auto pos1 = pos.begin();
		auto pos2 = pos.begin();
		advance(pos1, i-1);
		advance(pos2, i);
		Vector3d v1 = *pos1 - *pos0;
		Vector3d v2 = *pos2 - *pos0;
		Vector3d ret = v1.cross(v2);
		float v1L2 = v1.dot(v1);
		float v2L2 = v2.dot(v2);
		float retL2 = ret.dot(ret);
		if (v1L2 > 0 && v2L2 > 0 && retL2 / (v1L2 * v2L2) > EPS)
		{
			normalizeVec(ret);
			return ret;
		}
	}
	Vector3d ret(0,0,0);
	cout << "Not found face normal!"<<endl;
	return ret;
}

float getPolygonArea(list<Vector3d> pos)
{
	if (pos.size() < 3)
		return 0.0;
	auto pos0 = pos.begin();
	auto pos1 = pos.begin();
	advance(pos1, 1);
	Vector3d v1 = *pos1 - *pos0;
	Vector3d v2 = *pos1 - *pos0;
	float area = 0.0;
	for (int i = 2; i < pos.size(); ++i)
	{
		auto pos2 = pos.begin();
		v1 = v2;
		advance(pos2, i);
		v2 = *pos2 - *pos1;
		area = area + 0.5 * sqrt((v1.cross(v2).array() * v1.cross(v2).array()).sum());
	}
	return area;
}

Plant3D::Plant3D(Vector3d p0, Vector3d n)
{
	this->P0 = p0;
	this->N = n;
	normalizeVec(N);
	this->resetEquation();
}

void Plant3D::resetEquation()
{
	this->D = -P0.dot(N);
}

float Plant3D::disFromPlane(Vector3d P)
{
	return N.dot(P) + this->D;
}