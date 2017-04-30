#ifndef TRANSLATION 
#define TRANSLATION 

#include "geoTransform.h"

class translationCoords : public geoTransform {
public:
	float x, y, z;
	translationCoords(float fx, float fy, float fz) : x{fx}, y{fy}, z{fz} {};
	void apply() {
		glTranslatef(x, y, z);
	}

};

class translationTime : public geoTransform {
public:
	float time;
	std::vector<std::vector<float>> points;
	translationTime(float ftime) : time{ftime} {};

	void apply() {
		//glTranslatef(x, y, z);
	}

};
#endif
