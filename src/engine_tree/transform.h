#ifndef TRANSFORM
#define TRANSFORM

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

class transform {
public:
	virtual void apply() = 0;
};

#endif
