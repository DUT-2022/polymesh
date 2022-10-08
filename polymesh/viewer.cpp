#include <GL/glut.h>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <fstream>
#include "viewer.h"

using namespace std;
using namespace Eigen;

std::vector<Vector3d> points;
int faceNums;
float g_rotx = 0.0f, g_roty = 0.0f, g_rotz = 0.0f;
float g_modelPos[3] = { 0, 0, -15.0f };
float sx = 1.0f, sy = 1.0f, sz = 1.0f;
float g_scale = 10.0f;
const float g_lightPos[] = { 1.0f, 1.0f, 1.0f, 0.0f };


int readOff(std::string filepath, std::vector<Vector3d>& points)
{
	std::string line;
	std::ifstream fin(filepath);
	if (!fin.is_open())
	{
		std::cout << "文件 " << filepath << " 打开失败" << std::endl;
		exit(-1);
	}

	std::vector<Vector3d> vertexPosition; // off模型顶点位置和颜色
	fin >> line;    // OFF字符串
	// 读取顶点位置
	int vertexNum, faceNum, edgeNum;
	double x_max = -999999.0, x_min = 999999.0, y_max = -999999.0, y_min = 999999.0, z_max = -999999.0, z_min = 999999.0;
	fin >> vertexNum >> faceNum >> edgeNum;
	faceNums = faceNum;
	for (int i = 0; i < vertexNum; i++)
	{
		double p1, p2, p3;
		fin >> p1 >> p2 >> p3;
		vertexPosition.push_back(Vector3d(p1, p2, p3));
		if (p1 > x_max)
			x_max = p1;
		if (p1 < x_min)
			x_min = p1;
		if (p2 > y_max)
			y_max = p2;
		if (p2 < y_min)
			y_min = p2;
		if (p3 > z_max)
			z_max = p3;
		if (p3 < z_min)
			z_min = p3;
	}
	if (x_max - x_min > y_max - y_min)
	{
		if (x_max - x_min > z_max - z_min)
		{
			g_scale = g_scale / (x_max - x_min);
		}
		else
			g_scale = g_scale / (z_max - z_min);
	}
	else
	{
		if (y_max - y_min > z_max - z_min)
		{
			g_scale = g_scale / (y_max - y_min);
		}
		else
			g_scale = g_scale / (z_max - z_min);
	}
	for (int i = 0; i < vertexNum; i++)
	{
		vertexPosition[i][0] = g_scale * (vertexPosition[i][0] - 0.5 * (x_min + x_max));
		vertexPosition[i][1] = g_scale * (vertexPosition[i][1] - 0.5 * (y_min + y_max));
		vertexPosition[i][2] = g_scale * (vertexPosition[i][2] - 0.5 * (z_min + z_max));

	}
	// 根据面信息生成实际顶点
	points.clear();
	for (int i = 0; i < faceNum; i++)
	{
		int n, index1, index2, index3, index4;
		fin >> n >> index1 >> index2 >> index3 >> index4;
		points.push_back(vertexPosition[index1]);
		points.push_back(vertexPosition[index2]);
		points.push_back(vertexPosition[index3]);
		points.push_back(vertexPosition[index4]);
	}

}

int readMsh(std::string filepath, std::vector<Vector3d>& points)
{
	std::string line;
	std::ifstream fin(filepath);
	if (!fin.is_open())
	{
		std::cout << "文件 " << filepath << " 打开失败" << std::endl;
		exit(-1);
	}

	std::vector<Vector3d> vertexPosition; // off模型顶点位置和颜色
	fin >> line;    // OFF字符串
	// 读取顶点位置
	int vertexNum, faceNum, edgeNum;
	double x_max = -999999.0, x_min = 999999.0, y_max = -999999.0, y_min = 999999.0, z_max = -999999.0, z_min = 999999.0;
	
		

	vertexNum = 636;
	faceNum = 570;
	for (int i = 0; i < vertexNum; i++)
	{
		int p1;
		double p2, p3, p4;
		fin >> p1 >> p2 >> p3>> p4;
		cout << p1 << ' ' << p2 << ' ' << p3 << ' ' << p4 << endl;
		vertexPosition.push_back(Vector3d(p2, p3, p4));
		if (p1 > x_max)
			x_max = p1;
		if (p1 < x_min)
			x_min = p1;
		if (p2 > y_max)
			y_max = p2;
		if (p2 < y_min)
			y_min = p2;
		if (p3 > z_max)
			z_max = p3;
		if (p3 < z_min)
			z_min = p3;
	}
	if (x_max - x_min > y_max - y_min)
	{
		if (x_max - x_min > z_max - z_min)
		{
			g_scale = g_scale / (x_max - x_min);
		}
		else
			g_scale = g_scale / (z_max - z_min);
	}
	else
	{
		if (y_max - y_min > z_max - z_min)
		{
			g_scale = g_scale / (y_max - y_min);
		}
		else
			g_scale = g_scale / (z_max - z_min);
	}
	
	for (int i = 0; i < vertexNum; i++)
	{
		vertexPosition[i][0] = g_scale * (vertexPosition[i][0] - 0.5 * (x_min + x_max));
		vertexPosition[i][1] = g_scale * (vertexPosition[i][1] - 0.5 * (y_min + y_max));
		vertexPosition[i][2] = g_scale * (vertexPosition[i][2] - 0.5 * (z_min + z_max));

	}
	// 根据面信息生成实际顶点
	points.clear();

	faceNums = faceNum;
	for (int i = 0; i < faceNum; i++)
	{
		int n, index1, index2, index3, index4, index5, index6, index7, index8;
		fin >> n >> index1 >> index2 >> index3 >> index4 >> index5 >> index6 >> index7 >> index8;
		points.push_back(vertexPosition[index5]);
		points.push_back(vertexPosition[index6]);
		points.push_back(vertexPosition[index7]);
		points.push_back(vertexPosition[index8]);

	}

}

void DrawVertices()
{
	int i;

	glEnable(GL_LIGHTING);
	for (i = 0; i < faceNums; i++)
	{
		glBegin(GL_LINE_LOOP);
		{
			glVertex3d(points[4 * i + 0][0], points[4 * i][1], points[4 * i][2]);
			glVertex3d(points[4 * i + 1][0], points[4 * i + 1][1], points[4 * i + 1][2]);
			glVertex3d(points[4 * i + 2][0], points[4 * i + 2][1], points[4 * i + 2][2]);
			glVertex3d(points[4 * i + 3][0], points[4 * i + 3][1], points[4 * i + 3][2]);
		}
		glEnd();
	}
}

static void resize(int width, int height)
{
	const float ar = (float)width / (float)height;

	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45, ar, 1.0, 100.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

static void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glLoadIdentity();

	glTranslatef(g_modelPos[0], g_modelPos[1], -15);
	glRotatef(g_rotx, 1.0, 0.0, 0.0);
	glRotatef(g_roty, 0.0, 1.0, 0.0);
	glRotatef(g_rotz, 0.0, 0.0, 1.0);
	glScalef(sx, sy, 1);

	DrawVertices();

	glutSwapBuffers();
}


static void key(unsigned char key, int x, int y)
{
	int h;
	switch (key)
	{
	case 27:
		exit(0);
		break;
	case 'w':
		g_rotx -= 5;
		break;
	case 's':
		g_rotx += 5;
		break;
	case 'a':
		g_roty -= 5;
		break;
	case 'd':
		g_roty += 5;
		break;
	case 'q':
		g_rotz -= 5;
		break;
	case 'e':
		g_rotz += 5;
		break;
	case 'z':
		sx += 0.5;
		sy += 0.5;
		sz += 0.5;
		break;
	case 'x':
		sx -= 0.5;
		sy -= 0.5;
		sz -= 0.5;
		if (sx == 0)
			sx += 0.5;
		if (sy == 0)
			sy += 0.5;
		if (sz == 0)
			sz += 0.5;
		break;
	case '1':
		g_modelPos[0] += 0.1;
		break;
	case '3':
		g_modelPos[0] -= 0.1;
		break;
	case '2':
		g_modelPos[1] += 0.1;
		break;
	case '5':
		g_modelPos[1] -= 0.1;
		break;
	case 'o':
		g_rotx = 0.0f, g_roty = 0.0f, g_rotz = 0.0f;
		g_modelPos[0] = 0.0;
		g_modelPos[1] = 0.0;
		g_modelPos[2] = -15.0;
		sx = 1.0f, sy = 1.0f, sz = 1.0f;
		break;
	}

	glutPostRedisplay();
}

static void idle(void)
{
	glutPostRedisplay();
}

void initOff(int argc, char** argv, string filepath)
{
	glutInit(&argc, argv);
	glutInitWindowSize(640, 480);
	glutInitWindowPosition(10, 10);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);

	printf("Tips: use 'w', 's', 'a', 'd', 'z', 'x', 'q', 'e' keys to control the model\n");


	glutCreateWindow("off");

	glutReshapeFunc(resize);
	glutDisplayFunc(display);
	glutKeyboardFunc(key);
	glutIdleFunc(idle);

	glClearColor(1, 1, 1, 0);
	glDisable(GL_CULL_FACE);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);
	glLightfv(GL_LIGHT0, GL_POSITION, g_lightPos);

	glShadeModel(GL_SMOOTH);


	int faceNum = readOff(filepath, points);



	glutMainLoop();

}

