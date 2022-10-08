#include <iostream>
#include <string>
#include <stdio.h>
#include <io.h>
#include <cstring>
#include <fstream>
#include <math.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <gmsh.h_cwrap>
#include <list>

#include "CrossfieldCommon.h"
#include "PolyMesh.h"

using namespace std;
using namespace Eigen;

#define pi acos(-1)


int sign(double x)
{
	if (x > 0)
		return 1;
	else if (x == 0)
		return 0;
	else
		return -1;

}





Array3d arrayV(const Vector3d & v)
{
	//Transfer the input vector to an numpy array
	return Array3d(v(0), v(1), v(2));
}

MatrixXd matrixV(const Vector3d &v)
{
	//Transfter the input vector to a matrix
	MatrixXd xx(1, 3);
	xx << v(0), v(1), v(2);
	return xx;
}

MatrixXd matrixV_2(const Vector3d &v)
{
	//Transfter the input vector to a matrix

	return matrixV(v).transpose();

}

const Vector3d normaliseVec(const Vector3d &v)
{

	//Normalize a vector

	return v / v.norm();

}


list <MeshNode> flatten_list(list <list<MeshNode>> faceNodesList)
{
	list <MeshNode> FaceNodesList;
	for (auto List : faceNodesList)
	{
		for (auto node : List)
			FaceNodesList.push_back(node);
	}
	return FaceNodesList;
}


double dotProUnitVecs(const Vector3d &u1, const Vector3d &u2)
{
	//Dot product of input unit vectors
	double dotPro = u1.dot(u2);
	if (dotPro > 1)
		dotPro = 1;
	else if (dotPro < -1)
		dotPro = -1;

	return dotPro;

}

double extrapDistInTriElem(const Vector3d &p1, const Vector3d &p2, const Vector3d &p3,
	const double &d2, const double &d3, int localSize)
{
	//Calculate dis at p1 with d(x,y) as local linear approximation of Eikonal eq.
	//given d2, d3 | d1 > d2, d3
	const Vector3d v23 = p3 - p2;
	const Vector3d v21 = p1 - p2;

	double x3 = v23.norm();
	double l21 = v21.norm();

	double theta = acos(dotProUnitVecs(normaliseVec(v23), normaliseVec(v21)));
	double x1 = l21 * cos(theta);
	double y1 = l21 * sin(theta);

	if (y1 < 0)
		//raise IOError("y1<0");
		printf("IOError y1<0");
	double u3 = d3 - d2;

	double sqrtTerm = (1 - (u3 / x3)*(u3 / x3));
	if (sqrtTerm < 0)
	{
		printf("calcD: sqrtTerm<0");
		return NULL;
	}
	else
	{
		double u1 = (u3 / x3) * x1 + (sqrt(sqrtTerm))*y1;

		double d1 = d2 + u1 * localSize;
		return d1;
	}

}

bool isUnitVecBetweenUnitVecs(const Vector3d &u, const Vector3d &a, const Vector3d &b)
{
	//whether the input vector u is between two vector a and b

	//a and b have a subtended angle < pi
	if ((u.dot(a) == 1) || (u.dot(b) == 1))
		return true;
	double angab = acos(dotProUnitVecs(a, b));
	Vector3d crossProab = a.cross(b);
	double angau = acos(dotProUnitVecs(a, u));
	Vector3d crossProau = a.cross(u);

	if (sign(crossProab.dot(crossProau)) == -1)
		angau = (2 * 3.1415) - angau;
	double surplusAng = angab - angau;
	if (0 <= surplusAng && surplusAng <= angab)
		return true;
	else
		return false;

}

MatrixXd getA_R_B(const Vector3d &G_ia, const Vector3d &G_ja, const Vector3d &G_ka,
	const Vector3d &G_ib, const Vector3d &G_jb, const Vector3d &G_kb)
{
	//A rotation matrix
	MatrixXd A_R_B(3, 3);
	A_R_B << G_ia.dot(G_ib), G_ia.dot(G_jb), G_ia.dot(G_kb),
		G_ja.dot(G_ib), G_ja.dot(G_jb), G_ja.dot(G_kb),
		G_ka.dot(G_ib), G_ka.dot(G_jb), G_ka.dot(G_kb);

	return A_R_B;

}


double getAng(const Vector3d &cd1, const Vector3d &uref)
{
	//Get modulo 90 signed angle between 2 unit vectors
	double theta1 = (acos(dotProUnitVecs(uref, cd1))) / 3.1415 * 180;
	//theta = ((acos(dotProUnitVecs(uref, cd1)))/3.14*180) % 90.0;
	double theta = int(theta1) % 90;

	if (theta != 0 || theta != 180)
	{
		if (sign((uref.cross(cd1)).dot(Vector3d(0, 0, 1))) == -1)
			theta = -theta;
	}
	if (theta > 45)
		theta = -(90 - theta);

	else if (theta < -45)
		theta = (90 + theta);

	return theta;

}


Vector3d getClosestVec(const Vector3d &u1, const Matrix3d &uxs)
{

	//The the cloest vector in the uxs to u1
	double smallestAng = 1e9;
	Vector3d uxStar;
	int c = uxs.rows();
	for (int i = 0; i < c; i++)
	{
		Vector3d ux = uxs.row(i);
		double ang = acos(dotProUnitVecs(normaliseVec(u1), normaliseVec(ux)));
		if (ang < smallestAng)
		{
			smallestAng = ang;
			uxStar = ux;
		}
	}
	return uxStar;

}

Vector3d getPosIntersection2CoplanarTrajectories(const Vector3d &p1, const Vector3d &u1,
	const Vector3d &p2, const Vector3d &u2)
{
	double dotProu1u2 = dotProUnitVecs(u1, u2);
	if (dotProu1u2 == 1 || dotProu1u2 == -1 || (u1.cross(u2)).norm() == 0)
		//ruturn NULL;
		return Vector3d();

	if ((p2 - p1).norm() == 0)
		return p1;
	else
	{	//change basis
		//create local cartesian coord sys centred at p1
		//p1->q1, p2->q2
		Vector3d v12 = p2 - p1;
		Vector3d u12 = v12 / v12.norm();
		Vector3d ia = u12;
		Vector3d ka = u1.cross(u2) / (u1.cross(u2)).norm();
		Vector3d ja = ka.cross(ia);
		Matrix3d A_R_G = getA_R_B(ia, ja, ka, Vector3d(1, 0, 0), Vector3d(0, 1, 0), Vector3d(0, 0, 1));
		//G_R_A = A_R_G.transpose()

		Vector3d a_q1 = Vector3d(0, 0, 0);
		Vector3d a_q2 = Vector3d(A_R_G*matrixV(v12));
		Vector3d a_u1 = Vector3d(A_R_G*matrixV(u1));
		Vector3d a_u2 = Vector3d(A_R_G*matrixV(u2));

		double q1x = a_q1(0);
		double q1y = a_q1(1);
		double q2x = a_q2(0);
		double q2y = a_q2(1);
		double u1x = a_u1(0);
		double u1y = a_u1(1);
		double u2x = a_u2(0);
		double u2y = a_u2(1);

		double dotProu1u2 = dotProUnitVecs(Vector3d(u1x, u1y, 0), Vector3d(u2x, u2y, 0));
		if (dotProu1u2 == 1 || dotProu1u2 == -1 || (Vector3d((u1x, u1y, 0)).cross(Vector3d((u2x, u2y, 0)))).norm() == 0)
			//return NULL;
			return Vector3d();

		double t1 = (((q1x - q2x)*u2y) - ((q1y - q2y)*u2x)) / (u1y*u2x - u1x * u2y);
		Vector3d pint = p1 + t1 * u1;
		return pint;
	}
}

double getParaPosLine(const Vector3d &p0, const Vector3d &p1, const Vector3d &px)
{
	//Get the parameter of px in a line defined by p0, p1

	Vector3d v01 = p1 - p0;
	Vector3d v0x = px - p0;

	//project onto v01
	Vector3d u01 = v01 / v01.norm();
	Vector3d v0xProj = (v0x.dot(u01))*u01;

	double t;
	if (v0xProj.norm() == 0)
		t = 0;
	else
		t = (v0xProj.norm() / v01.norm())*(sign(v0xProj.dot(v01)));

	return t;

}

Vector2d lineIntersection(const Vector3d &p1, const Vector3d &p2, const Vector3d &pa, const Vector3d &pb)
{
	//Calculate the intersection of two 2D lines, defined by p1 p2 and pa pb

	Vector3d L1 = Vector3d(p1(1) - p2(1), p2(0) - p1(0), p1(0) * p2(1) - p2(0) * p1(1));
	Vector3d L2 = Vector3d(pa(1) - pb(1), pb(0) - pa(0), pa(0) * pb(1) - pb(0) * pa(1));

	double D = L1(0) * L2(1) - L1(1) * L2(0);
	double Dx = L1(2) * L2(1) - L1(1) * L2(2);
	double Dy = L1(0) * L2(2) - L1(2) * L2(0);
	if (D != 0)
	{
		double x = Dx / D;
		double y = Dy / D;
		return Vector2d(x, y);
	}
	else
		//return NULL;
		return Array2d();
}

double HertzmannZorinDistortion(double theta1, double theta2)
{
	//HertzmannZorin Distortion, see thesis

	return 0.5*(1.0 - cos(4 * theta1 - 4 * theta2));

}

bool isPntOnLineSeg(const Vector3d &a, const Vector3d &b, const Vector3d &c)
{
	//Determine if a point c lies between point a and point b(already know c is on the
	//line passing a and b)

	double dotproduct = (c(0) - a(0)) * (b(0) - a(0)) + (c(1) - a(1))*(b(1) - a(1));
	if (dotproduct < 0)
		return false;

	double squaredlengthba = (b(0) - a(0))*(b(0) - a(0)) + (b(1) - a(1))*(b(1) - a(1));
	if (dotproduct > squaredlengthba)
		return false;

	return true;

}


Vector3d projVecOntoPlane(const Vector3d &v, const Vector3d &n)
{
	//Project a vector onto a plane

	Vector3d u = normaliseVec(v);
	if ((n.cross(u)).norm() < 0.0001)
		//raise IOError("norm(cross(n,u))<0.0001");
		printf("IOError  norm(cross(n,u))<0.0001");

	Vector3d kq = n;
	Vector3d jq = normaliseVec(kq.cross(u));
	Vector3d iq = jq.cross(kq);

	Matrix3d G_R_Q = getA_R_B(Vector3d(1, 0, 0), Vector3d(0, 1, 0), Vector3d(0, 0, 1), iq, jq, kq);
	Matrix3d Q_R_G = G_R_Q.transpose();
	Vector3d q_v = Vector3d(Q_R_G*matrixV_2(v));
	Vector3d q_vnew = Vector3d(q_v(0), q_v(1), 0);
	Vector3d vnew = Vector3d(G_R_Q*matrixV_2(q_vnew));
	return vnew;

}


MatrixXd RotMatrix(const Vector3d &n, const double &angle)   //angle..degrees
{
	//Rotating a vector v around a vector n by the input angle to get a new vector VV.
	//The new vector VV can be calculated by VV = R * v, where R is the rotation metrix
	//calculated in this function.
	//https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
	Matrix3d R;
	if (angle == 0)
		R << 1, 0, 0, 0, 1, 0, 0, 0, 1;
	else
	{
		double theta = angle * (3.14 / 180);
		double R11 = cos(theta) + (n(0) * n(0))*(1 - cos(theta));
		double R12 = n(0) * n(1) * (1 - cos(theta)) - n(2) * sin(theta);
		double R13 = n(0) * n(2) * (1 - cos(theta)) + n(1) * sin(theta);
		double R21 = n(1) * n(0) * (1 - cos(theta)) + n(2) * sin(theta);
		double R22 = cos(theta) + (n(1) * n(1))*(1 - cos(theta));
		double R23 = n(1) * n(2) * (1 - cos(theta)) - n(0) * sin(theta);
		double R31 = n(2) * n(0) * (1 - cos(theta)) - n(1) * sin(theta);
		double R32 = n(2) * n(1) * (1 - cos(theta)) + n(0) * sin(theta);
		double R33 = cos(theta) + (n(2) * n(2))*(1 - cos(theta));

		R << R11, R12, R13, R21, R22, R23, R31, R32, R33;
	}
	return R;

}


Vector3d rotVecFromPlaneBToPlaneA(const Vector3d &v, const Vector3d &na, const Vector3d &nb)
{
	//Rotate a vector from plane B to plane A

	Vector3d crossnbna = nb.cross(na);
	Vector3d vrot;

	if ((crossnbna).norm() > 0.001)
	{
		double theta = (acos(dotProUnitVecs(nb, na))) / 3.14 * 180;
		Vector3d u = normaliseVec(crossnbna);
		Matrix3d RotMat = RotMatrix(u, theta);
		vrot = normaliseVec(arrayV(RotMat*matrixV_2(v)));
	}
	else
	{//RotMat = matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
		vrot = v;
	}
	return vrot;

}


MatrixXd RotMatrix2(const Vector3d &rv)    //..radians
{
	Matrix3d R;
	double theta = rv.norm();
	if (theta == 0)
		R << 1, 0, 0, 0, 1, 0, 0, 0, 1;
	else
	{
		Vector3d n = rv / rv.norm();
		double R11 = cos(theta) + (n(0) * n(0))*(1 - cos(theta));
		double R12 = n(0) * n(1) * (1 - cos(theta)) - n(2) * sin(theta);
		double R13 = n(0) * n(2) * (1 - cos(theta)) + n(1) * sin(theta);
		double R21 = n(1) * n(0) * (1 - cos(theta)) + n(2) * sin(theta);
		double R22 = cos(theta) + (n(1) * n(1))*(1 - cos(theta));
		double R23 = n(1) * n(2) * (1 - cos(theta)) - n(0) * sin(theta);
		double R31 = n(2) * n(0) * (1 - cos(theta)) - n(1) * sin(theta);
		double R32 = n(2) * n(1) * (1 - cos(theta)) + n(0) * sin(theta);
		double R33 = cos(theta) + (n(2) * n(2))*(1 - cos(theta));
		R << R11, R12, R13, R21, R22, R23, R31, R32, R33;
	}
	return R;

}


MatrixXd getEllipseVectorsOnPlane(const Matrix3d &tensor, const Vector3d &planeNormal)
{
	//get coordinate frame with k aligned with n
	Vector3d kq = planeNormal;
	Vector3d iq;
	//print "kq", kq
	if (abs(kq.dot(Vector3d(1, 0, 0))) < 0.99)
		iq = normaliseVec(kq.cross(Vector3d(1, 0, 0)));
	else
		iq = normaliseVec(kq.cross(Vector3d(0, 1, 0)));
	Vector3d jq = kq.cross(iq);

	Matrix3d G_R_Q = getA_R_B(Vector3d(1, 0, 0), Vector3d(0, 1, 0), Vector3d(0, 0, 1), iq, jq, kq);
	Matrix3d Q_R_G = G_R_Q.transpose();

	Matrix3d Adash = Q_R_G * tensor* G_R_Q;
	Matrix2d Alpha;
	Alpha << Adash(0, 0), Adash(0, 1), Adash(0, 1), Adash(1, 1);

	EigenSolver<MatrixXd> es(Alpha);
	MatrixXd evalues = es.pseudoEigenvalueMatrix();
	MatrixXd EVECS = es.pseudoEigenvectors();

	double l1 = evalues(0, 0);
	double l2 = evalues(1, 1);
	Vector2d evec1 = EVECS.col(0);
	Vector2d evec2 = EVECS.col(1);
	Vector2d a1 = evec1 * (pow(l1, -0.5)) / 2.0;
	Vector2d a2 = evec2 * (pow(l2, -0.5)) / 2.0;
	Vector3d q_a1 = Vector3d(a1(0), a1(1), 0);
	Vector3d q_a2 = Vector3d(a2(0), a2(1), 0);
	Vector3d g_a1 = Vector3d(G_R_Q*matrixV_2(q_a1));
	Vector3d g_a2 = Vector3d(G_R_Q*matrixV_2(q_a2));
	MatrixXd g(2, 3);
	g << g_a1(0), g_a1(1), g_a1(2), g_a2(0), g_a2(1), g_a2(2);
	return g;

}


Vector3d lineIntersetOnPlane(const Vector3d &p1, const Vector3d &u1,
	const Vector3d &p2, const Vector3d &u2)
{
	//Calculate the intersection btwwen a line defined by point p1 and a vector u1
	//and another line defined by p2 and a vector u2.u1 and u2 are on the same plane.

	double dotProu1u2 = dotProUnitVecs(u1, u2);
	if (dotProu1u2 == 1 || dotProu1u2 == -1 || (u1.cross(u2)).norm() == 0)
		//return NULL;
		return Vector3d();

	if ((p2 - p1).norm() == 0)
		return p1;
	else
	{
		//change basis
		//create local cartesian coord sys centred at p1
		//p1->q1, p2->q2
		Vector3d v12 = p2 - p1;
		Vector3d u12 = v12 / v12.norm();
		Vector3d ia = u12;
		Vector3d ka = u1.cross(u2) / (u1.cross(u2)).norm();
		Vector3d ja = ka.cross(ia);
		Matrix3d A_R_G = getA_R_B(ia, ja, ka, Vector3d(1, 0, 0), Vector3d(0, 1, 0), Vector3d(0, 0, 1));
		//G_R_A = A_R_G.transpose()

		Vector3d a_q1 = Vector3d(0, 0, 0);
		Vector3d a_q2 = Vector3d(A_R_G*matrixV_2(v12));
		Vector3d a_u1 = Vector3d(A_R_G*matrixV_2(u1));
		Vector3d a_u2 = Vector3d(A_R_G*matrixV_2(u2));

		double q1x = a_q1(0);
		double q1y = a_q1(1);
		double q2x = a_q2(0);
		double q2y = a_q2(1);
		double u1x = a_u1(0);
		double u1y = a_u1(1);
		double u2x = a_u2(0);
		double u2y = a_u2(1);

		double dotProu1u2 = dotProUnitVecs(Vector3d(u1x, u1y, 0), Vector3d(u2x, u2y, 0));
		if (dotProu1u2 == 1 || dotProu1u2 == -1 || ((Vector3d(u1x, u1y, 0)).cross(Vector3d(u2x, u2y, 0))).norm() == 0)
			//return NULL;
			return Vector3d();

		double	t1 = (((q1x - q2x)*u2y) - ((q1y - q2y)*u2x)) / (u1y*u2x - u1x * u2y);
		Vector3d pint = p1 + t1 * u1;
		return pint;

	}
}