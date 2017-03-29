#ifndef SCALE
#define SCALE

#include "transform.h"

class scale : public transform {
public:
	float x, y, z;
	scale(float fx, float fy, float fz) : x{fx}, y{fy}, z{fz} {};
	void apply() {
		glScalef(x, y, z);
	}

};
#endif
