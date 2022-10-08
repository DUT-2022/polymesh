#pragma once
#include <iostream>
#include <math.h>
#define PI acos(-1)

using namespace std;
class CrossFieldParams
{
public:
	//Background tri mesh size
	int bgMeshSize = 8;

	//Scaling factor of coefficient w1
	int w1Scale = 1000000;

	//Scaling factor of coefficient w2
	int w2Scale = 100;

	//Smoothing iteration number
	int nIter = 5;

	//
	double maxSmooAng = 15 / 180 * PI;

	//To avoid a pair of positive and negative singulariteis
	//appear in one background element
	bool avoidDipolesBool = false;

	//Generate partial blocks insteand of full blocks
	bool partialBlockBool = false;
};
