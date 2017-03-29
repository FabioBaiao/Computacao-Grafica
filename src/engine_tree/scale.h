#ifndef SCALE
#define SCALE

#include "geoTransform.h"

class scale : public geoTransform {
public:
	float x, y, z;
	scale(float fx, float fy, float fz) : x{fx}, y{fy}, z{fz} {};
	void apply() {
		glScalef(x, y, z);
	}

};
#endif
