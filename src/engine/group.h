#ifndef GROUP
#define GROUP

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "geoTransform.h"

class color {
public:
	float colors[4];
	GLenum component;
	color(float rr, float gg, float bb, GLenum ecomponent) : 
		colors{rr, gg, bb, 1.0f}, component{ecomponent} {};

	/*void apply(){
		glMaterialfv(GL_FRONT, component, colors);
	}*/
};

class randSpecs {
public:
	int n;
	float minR, maxR, minS, maxS;
	randSpecs(int nn, float minRR, float maxRR, float minSS, float maxSS) : 
		n{nn}, minR{minRR}, maxR{maxRR}, minS{minSS}, maxS{maxSS} {};
};

class randomModel {
public:
	std::vector<std::string> models;
	std::vector<std::vector<color>> modelsColor;
	std::vector<randSpecs> specs;
};

class group {
public:
	std::vector<std::string> models;
	std::vector<std::vector<color>> modelsColor;
	std::vector<GLuint> modelsTextures;
	std::vector<geoTransform *> transforms; 
	std::vector<group> childGroups;
	std::vector<randomModel> randoms;
};

#endif
