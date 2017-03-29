#ifndef ROTATE 
#define ROTATE 

#include "transform.h"

class rotate : public transform {
public:
	float angle, x, y, z;
	rotate(float fangle, float fx, float fy, float fz) : angle{fangle}, x{fx}, y{fy}, z{fz} {};
	
	void apply() {
		glRotatef(angle, x, y, z);
	}

};
#endif
