#ifndef POINT
#define POINT

// Imutable three-dimensional points
class point{
public:
	const float x, y, z;
	point(float px, float py, float pz):x{px}, y{py}, z{pz} {};
	point():x{}, y{}, z{} {};
	~point(){};
};

#endif //ifndef POINT
