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
	float referential[16];
	std::vector<figure> figures;
};

#endif
