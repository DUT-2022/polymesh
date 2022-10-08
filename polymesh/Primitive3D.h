#pragma once
#include <iostream>
#include <list>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

auto normalizeVec(Vector3d&);
Vector3d getFaceNormal(list<Vector3d>);
float getPolygonArea(list<Vector3d>);
class Plant3D
{
public:
	Vector3d P0;
	Vector3d N;
	float D = 0;
	Plant3D(Vector3d, Vector3d);
	void resetEquation();
	float disFromPlane(Vector3d);
};

