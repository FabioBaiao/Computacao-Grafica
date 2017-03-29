#ifndef GROUP
#define GROUP

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "point.h"
#include "triangle.h"
#include "transform.h"

class color{
public:
	float r, g, b;
	color(float rr, float gg, float bb) : r{rr}, g{gg}, b{bb} {};
};

class group{
public:
	std::vector<std::string> models;
	std::vector<color> models_color;
	std::vector<transform*> transforms; 
	std::vector<group> child_groups;
};

#endif
