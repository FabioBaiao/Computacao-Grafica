#ifndef GROUP
#define GROUP

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "point.h"
#include "triangle.h"


class group{
public:
	//string model_name;
	//std::vector<figure> figures;
	std::vector<string> models;
	//float referential[16];
	std::vector<transform> transforms; // fazer class ou por já campos das transformacoes? 
	group * child_group;
};

#endif
