#ifndef ROTATION 
#define ROTATION 

#include "geoTransform.h"

class rotationAngle : public geoTransform {
public:
	float angle, x, y, z;
	rotationAngle(float fangle, float fx, float fy, float fz) : angle{fangle}, x{fx}, y{fy}, z{fz} {};
	
	void apply() {
		glRotatef(angle, x, y, z);
	}

};

class rotationTime : public geoTransform {
public:
	float time, x, y, z;
	rotationTime(float ftime, float fx, float fy, float fz) : time{ftime}, x{fx}, y{fy}, z{fz} {};
	
	void apply() {
		//glRotatef(time, x, y, z);
	}

};
#endif
