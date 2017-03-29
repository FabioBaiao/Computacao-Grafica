#ifndef TRANSLATE 
#define TRANSLATE 

#include "transform.h"

class translate : public transform {
public:
	float x, y, z;
	translate(float fx, float fy, float fz) : x{fx}, y{fy}, z{fz} {};
	void apply() {
		glTranslatef(x, y, z);
	}

};
#endif
