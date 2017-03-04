//SORRY FOR BAD ENGLISH :/

/*
	The final version are the three drawPlane__ functions
	-----------------------------------------------------
*/

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
	if 'orient' is 1, the plane will be oriented upwards. if it is -1, will be oriented downwards
*/

std::string drawPlaneXZ (float x, float y, float z, int nDiv, int orient){
	
	std::string ret = "";

	int count = 0;

	float deltaX = x/nDiv;
	float deltaZ = z/nDiv;

	x = x/2;
	z = z/2;

	float xi = x;
	float zi = z;

	for (int i = 1; i <= nDiv; i++){
		for (int j = 1; j <= nDiv; j++){

			ret += std::to_string(xi) + " " + std::to_string(y) + " " + std::to_string(zi) + "\n";
			if (orient == 1)
				ret += std::to_string(xi) + " " + std::to_string(y) + " " + std::to_string(zi-deltaZ) + "\n";
			else if (orient == -1)
				ret += std::to_string(xi-deltaX) + " " + std::to_string(y) + " " + std::to_string(zi) + "\n";

			ret += std::to_string(xi-deltaX) + " " + std::to_string(y) + " " + std::to_string(zi-deltaZ) + "\n";

			ret += std::to_string(xi) + " " + std::to_string(y) + " " + std::to_string(zi) + "\n";
			ret += std::to_string(xi-deltaX) + " " + std::to_string(y) + " " + std::to_string(zi-deltaZ) + "\n";
			if (orient == 1)
				ret += std::to_string(xi-deltaX) + " " + std::to_string(y) + " " + std::to_string(zi) + "\n";
			else if (orient == -1)
				ret += std::to_string(xi) + " " + std::to_string(y) + " " + std::to_string(zi-deltaZ) + "\n";

			count += 6;
			xi = x - j * deltaX;
		}

		zi = z - i * deltaZ;
		xi = x;
	}
	return std::to_string(count) + "\n" + ret;
}

std::string drawPlaneXY (float x, float y, float z, int nDiv, int orient){
	
	std::string ret = "";

	int count = 0;

	float deltaX = x/nDiv;
	float deltaY = y/nDiv;

	x = x/2;
	y = y/2;

	float xi = x;
	float yi = -y;

	for (int i = 1; i <= nDiv; i++){
		for (int j = 1; j <= nDiv; j++){

			ret += std::to_string(xi) + " " + std::to_string(yi) + " " + std::to_string(z) + "\n";
			if (orient == 1)
				ret += std::to_string(xi) + " " + std::to_string(yi+deltaY) + " " + std::to_string(z) + "\n";
			else if (orient == -1)
				ret += std::to_string(xi+deltaX) + " " + std::to_string(yi) + " " + std::to_string(z) + "\n";

			ret += std::to_string(xi-deltaX) + " " + std::to_string(yi+deltaY) + " " + std::to_string(z) + "\n";

			ret += std::to_string(xi) + " " + std::to_string(yi) + " " + std::to_string(z) + "\n";
			ret += std::to_string(xi-deltaX) + " " + std::to_string(yi+deltaY) + " " + std::to_string(z) + "\n";
			if (orient == 1)
				ret += std::to_string(xi-deltaX) + " " + std::to_string(yi) + " " + std::to_string(z) + "\n";
			else if (orient == -1)
				ret += std::to_string(xi) + " " + std::to_string(yi-deltaY) + " " + std::to_string(z) + "\n";

			count += 6;
			xi = x - j * deltaX;
		}

		yi = -y + i * deltaY;
		xi = x;
	}
	return std::to_string(count) + "\n" + ret;
}

std::string drawPlaneYZ (float x, float y, float z, int nDiv, int orient){
	
	std::string ret = "";

	int count = 0;

	float deltaY = y/nDiv;
	float deltaZ = z/nDiv;

	y = y/2;
	z = z/2;

	float yi = -y;
	float zi = -z;

	for (int i = 1; i <= nDiv; i++){
		for (int j = 1; j <= nDiv; j++){

			ret += std::to_string(x) + " " + std::to_string(yi) + " " + std::to_string(zi) + "\n";
			if (orient == 1)
				ret += std::to_string(x) + " " + std::to_string(yi+deltaY) + " " + std::to_string(zi) + "\n";
			else if (orient == -1)
				ret += std::to_string(x) + " " + std::to_string(yi) + " " + std::to_string(zi+deltaZ) + "\n";

			ret += std::to_string(x) + " " + std::to_string(yi+deltaY) + " " + std::to_string(zi+deltaZ) + "\n";

			ret += std::to_string(x) + " " + std::to_string(yi) + " " + std::to_string(zi) + "\n";
			ret += std::to_string(x) + " " + std::to_string(yi+deltaY) + " " + std::to_string(zi+deltaZ) + "\n";
			if (orient == 1)
				ret += std::to_string(x) + " " + std::to_string(yi) + " " + std::to_string(zi+deltaZ) + "\n";
			else if (orient == -1)
				ret += std::to_string(x) + " " + std::to_string(yi+deltaY) + " " + std::to_string(zi) + "\n";

			count += 6;
			zi = z + j * deltaZ;
		}

		yi = y + i * deltaY;
		zi = z;
	}
	return std::to_string(count) + "\n" + ret;
}

/*RIGHT HAND RULE

	X -> 90º rotation
	2X -> 180º rotation
	3X -> 270º/-90º rotation

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
*/