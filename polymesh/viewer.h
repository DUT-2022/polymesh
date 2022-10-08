#pragma once
#include <GL/glut.h>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <fstream>

using namespace std;
using namespace Eigen;

int readOff(std::string, std::vector<Vector3d>&);
int readMsh(std::string, std::vector<Vector3d>&);
void DrawVertices();
static void resize(int , int );
static void display(void);
static void key(unsigned char, int, int);
static void idle(void);
void initOff(int , char** , string);
