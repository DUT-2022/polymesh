#pragma once
#include <iostream>
#include "PolyMesh.h"
using namespace std;

int sign(double x);
Array3d arrayV(const Vector3d & v);
MatrixXd matrixV(const Vector3d &v);
MatrixXd matrixV_2(const Vector3d &v);
const Vector3d normaliseVec(const Vector3d &v);
list <MeshNode> flatten_list(list <list<MeshNode>> faceNodesList);
double dotProUnitVecs(const Vector3d &u1, const Vector3d &u2);
double extrapDistInTriElem(const Vector3d &p1, const Vector3d &p2, const Vector3d &p3,
	const double &d2, const double &d3, int localSize = 1);
bool isUnitVecBetweenUnitVecs(const Vector3d &u, const Vector3d &a, const Vector3d &b);
MatrixXd getA_R_B(const Vector3d &G_ia, const Vector3d &G_ja, const Vector3d &G_ka,
	const Vector3d &G_ib, const Vector3d &G_jb, const Vector3d &G_kb);
double getAng(const Vector3d &cd1, const Vector3d &uref);
Vector3d getClosestVec(const Vector3d &u1, const Matrix3d &uxs);
Vector3d getPosIntersection2CoplanarTrajectories(const Vector3d &p1, const Vector3d &u1,
	const Vector3d &p2, const Vector3d &u2);
double getParaPosLine(const Vector3d &p0, const Vector3d &p1, const Vector3d &px);
Vector2d lineIntersection(const Vector3d &p1, const Vector3d &p2, const Vector3d &pa, const Vector3d &pb);
double HertzmannZorinDistortion(double theta1, double theta2);
bool isPntOnLineSeg(const Vector3d &a, const Vector3d &b, const Vector3d &c);
Vector3d projVecOntoPlane(const Vector3d &v, const Vector3d &n);
MatrixXd RotMatrix(const Vector3d &n, const double &angle);   //angle..degrees
Vector3d rotVecFromPlaneBToPlaneA(const Vector3d &v, const Vector3d &na, const Vector3d &nb);
MatrixXd RotMatrix2(const Vector3d &rv);    //..radians;
MatrixXd getEllipseVectorsOnPlane(const Matrix3d &tensor, const Vector3d &planeNormal);
Vector3d lineIntersetOnPlane(const Vector3d &p1, const Vector3d &u1,
	const Vector3d &p2, const Vector3d &u2);
