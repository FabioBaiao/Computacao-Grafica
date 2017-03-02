//SORRY FOR BAD ENGLISH :/

/*
	This is for the simplest case of the plane
	When doing the other cases, they can be used for this case and this one can be removed
*/

void plane (float l, FILE* f){

	/* 
		l must be greater than 0:
		validate here or before?
	*/

	l = l/2;

	/*
		print the following to f

	vertices:
	A -> l 0 l
	B -> l 0 -l
	C -> -l 0 -l
	D -> -l 0 l
	//number of vertices for the triangles (2 triangles each with 3 vertices)
	6
	//triangle #1
	vertex A
	vertex B
	vertex D
	//triangle #2
	vertex B
	vertex C
	vertex D
	*/
}

/*
	Calculates the points to draw a plane paralel to the x0z face.
	The 'y' argument is to allow to draw a plane out of the x0z face.
*/

std::string plane2 (float x, float y, float z, int nDiv){
	
	std::string ret = "";

	int count = 0;

	float deltaX = x/nDiv;
	float deltaZ = z/nDiv;

	x = x/2;
	z = z/2;

	for (float zi = z; zi > -z; zi -= deltaZ){
		for (float xi = x; x > -x; xi -= deltaX){

			ret += std::to_string(xi) + " " + std::to_string(y) + " " + std::to_string(zi) + "\n";
			ret += std::to_string(xi) + " " + std::to_string(y) + " " + std::to_string(zi-deltaZ) + "\n";
			ret += std::to_string(xi-deltaX) + " " + std::to_string(y) + " " + std::to_string(zi-deltaZ) + "\n";

			ret += std::to_string(xi) + " " + std::to_string(y) + " " + std::to_string(zi) + "\n";
			ret += std::to_string(xi-deltaX) + " " + std::to_string(y) + " " + std::to_string(zi-deltaZ) + "\n";
			ret += std::to_string(xi-deltaX) + " " + std::to_string(y) + " " + std::to_string(zi) + "\n";

			count += 6;
		}
	}
	return std::to_string(count) + "\n" + ret;
}

/*RIGHT HAND RULE

	X -> 90º rotation
	2X -> 180º rotation
	3X -> 270º/-90º rotation
*/

std::string rotateX (std::string vertex){
	float x, y, z;
	x = ...;
	y = ...;
	z = ...;
	float aux = y;
	y = -z;
	z = aux;
	return std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + "\n";
}

std::string rotate2X (std::string vertex){
	float x, y, z;
	x = ...;
	y = ...;
	z = ...;
	y = -y;
	z = -z;
	return std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + "\n";
}

std::string rotate3X (std::string vertex){
	float x, y, z;
	x = ...;
	y = ...;
	z = ...;
	float aux = y;
	y = z;
	z = -aux;
	return std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + "\n";
}

std::string rotateZ (std::string vertex){
	float x, y, z;
	x = ...;
	y = ...;
	z = ...;
	float aux = x;
	x = -y;
	y = aux;
	return std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + "\n";
}

std::string rotate3Z (std::string vertex){
	float x, y, z;
	x = ...;
	y = ...;
	z = ...;
	float aux = x;
	x = y;
	y = -aux;
	return std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + "\n";
}