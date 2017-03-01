#ifndef POINT
#define POINT

// Imutable three-dimensional points
class point{
	const float x, y, z;
public:
	point(float px, float py, float pz):x{px}, y{py}, z{pz} {};
};

#endif //ifndef POINT
