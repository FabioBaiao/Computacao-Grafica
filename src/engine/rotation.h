#ifndef ROTATION 
#define ROTATION 

#include "geoTransform.h"

class rotation : public geoTransform {
public:
	float angle, x, y, z;
	rotation(float fangle, float fx, float fy, float fz) : angle{fangle}, x{fx}, y{fy}, z{fz} {};
	
	void apply() {
		glRotatef(angle, x, y, z);
	}

};
#endif
