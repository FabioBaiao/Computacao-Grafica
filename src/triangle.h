#ifndef TRIANGLE
#define TRIANGLE
#include "point.h"

class triangle{
public:
	point p1, p2, p3;
	float color_r, color_g, color_b;
	triangle(point pt1, point pt2, point pt3, float c_r, float c_g,
		 float c_b):p1{pt1}, p2{pt2}, p3{pt3}, color_r{c_r}, color_g{c_g}, color_b{c_b} {};
	triangle(point pt1, point pt2, point pt3):p1{pt1}, p2{pt2}, p3{pt3}, color_r{0}, color_g{0}, color_b{0} {};
};

#endif
