#ifndef TRANSLATION 
#define TRANSLATION 

#include "geoTransform.h"

class translation : public geoTransform {
public:
	float x, y, z;
	translation(float fx, float fy, float fz) : x{fx}, y{fy}, z{fz} {};
	void apply() {
		glTranslatef(x, y, z);
	}

};
#endif
