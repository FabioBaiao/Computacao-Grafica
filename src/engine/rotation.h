#ifndef ROTATION 
#define ROTATION 

#include "geoTransform.h"

class rotation : public geoTransform {
public:
	float angle, x, y, z, time;
	rotation(float fangle, float fx, float fy, float fz, float ftime) : angle{fangle}, x{fx}, y{fy}, z{fz}, time{ftime*(float)1e3} {};
	
	void apply() {
		angle = (time == 0 ? angle : glutGet(GLUT_ELAPSED_TIME)*360/time);
		glRotatef(angle, x, y, z);
	}

};
#endif
