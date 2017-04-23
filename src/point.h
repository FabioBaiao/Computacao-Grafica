#ifndef POINT
#define POINT

// Three-dimensional point
class point{
public:
	float x, y, z;
	point(float px, float py, float pz):x{px}, y{py}, z{pz} {};
	point(const point &p):x{p.x}, y{p.y}, z{p.z} {};
	point():x{}, y{}, z{} {};
	~point(){};
};

#endif //ifndef POINT
