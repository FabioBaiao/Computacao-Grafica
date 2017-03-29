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
	std::vector<figure> figures;
	float referential[16];
};

#endif
