#ifndef GEO_TRANSFORM
#define GEO_TRANSFORM

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

class geoTransform {
public:
	virtual void apply() = 0;
};

#endif
