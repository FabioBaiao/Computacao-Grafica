std::string planeXZ (float x, float y, float z, int nDiv, int orient){
	
	std::ostringstream os;

	int count = 0;

	float deltaX = x/nDiv;
	float deltaZ = z/nDiv;

	x = x/2;
	z = z/2;

	float xi = x;
	float zi = z;

	for (int i = 1; i <= nDiv; i++){
		for (int j = 1; j <= nDiv; j++){

			if (orient == 1){
				os << xi << " " << y << " " << zi << endl;				
				os << xi << " " << y << " " << (zi-deltaZ) << endl;
				os << (xi-deltaX) << " " << y << " " << (zi-deltaZ) << endl;

				os << xi << " " << y << " " << zi << endl;				
				os << (xi-deltaX) << " " << y << " " << (zi-deltaZ) << endl;
				os << (xi-deltaX) << " " << y << " " << zi << endl;
			}
			else {
				os << xi << " " << y << " " << zi << endl;				
				os << (xi-deltaX) << " " << y << " " << zi << endl;
				os << (xi-deltaX) << " " << y << " " << (zi-deltaZ) << endl;

				os << xi << " " << y << " " << zi << endl;				
				os << (xi-deltaX) << " " << y << " " << (zi-deltaZ) << endl;
				os << xi << " " << y << " " << (zi-deltaZ) << endl;
			}

			count += 6;
			xi = x - j * deltaX;
		}

		zi = z - i * deltaZ;
		xi = x;
	}
	return std::to_string(count) + "\n" + os.str();
}

std::string planeXY (float x, float y, float z, int nDiv, int orient){
	
	std::ostringstream os;

	int count = 0;

	float deltaX = x/nDiv;
	float deltaY = y/nDiv;

	x = x/2;
	y = y/2;

	float xi = x;
	float yi = -y;

	for (int i = 1; i <= nDiv; i++){
		for (int j = 1; j <= nDiv; j++){

			if (orient == 1){
				os << xi << " " << yi << " " << z << endl;
				os << xi << " " << (yi+deltaY) << " " << z << endl;
				os << (xi-deltaX) << " " << (yi+deltaY) << " " << z << endl; 

				os << xi << " " << yi << " " << z << endl;
				os << (xi-deltaX) << " " << (yi+deltaY) << " " << z << endl; 				
				os << (xi-deltaX) << " " << yi << " " << z << endl;
			}
			else {
				os << xi << " " << yi << " " << z << endl;
				os << (xi-deltaX) << " " << yi << " " << z << endl;	
				os << (xi-deltaX) << " " << (yi+deltaY) << " " << z << endl; 

				os << xi << " " << yi << " " << z << endl;
				os << (xi-deltaX) << " " << (yi+deltaY) << " " << z << endl; 				
				os << xi << " " << (yi+deltaY) << " " << z << endl;
			}

			count += 6;
			xi = x - j * deltaX;
		}

		yi = -y + i * deltaY;
		xi = x;
	}
	return std::to_string(count) + "\n" + os.str();
}

std::string planeYZ (float x, float y, float z, int nDiv, int orient){
	
	std::ostringstream;

	int count = 0;

	float deltaY = y/nDiv;
	float deltaZ = z/nDiv;

	y = y/2;
	z = z/2;

	float yi = -y;
	float zi = -z;

	for (int i = 1; i <= nDiv; i++){
		for (int j = 1; j <= nDiv; j++){

			if (orient == 1){
				os << x << " " << yi << " " << zi << endl;
				os << x << " " << (yi+deltaY) << " " << zi << endl;
				os << x << " " << (yi+deltaY) << " " << (zi+deltaZ) << endl;

				os << x << " " << yi << " " << zi << endl;
				os << x << " " << (yi+deltaY) << " " << (zi+deltaZ) << endl;
				os << x << " " << yi << " " << (zi+deltaZ) << endl;
			}
			else {
				os << x << " " << yi << " " << zi << endl;
				os << x << " " << yi << " " << (zi+deltaZ) << endl;
				os << x << " " << (yi+deltaY) << " " << (zi+deltaZ) << endl;
				
				os << x << " " << yi << " " << zi << endl;
				os << x << " " << (yi+deltaY) << " " << (zi+deltaZ) << endl;
				os << x << " " << (yi+deltaY) << " " << zi << endl;

			}

			count += 6;
			zi = -z + j * deltaZ;
		}

		yi = -y + i * deltaY;
		zi = -z;
	}
	return std::to_string(count) + "\n" + os.str();
}
