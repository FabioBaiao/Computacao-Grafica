std::string box (float x, float y, float z, int nDiv){
	std::string plane, line;
	std::ostringstream os;
	std::istringstream iss;
	int vertices, count;

	plane = drawPlaneXY(x, y, z/2, nDiv, 1);
	std::istringstream iss1 (plane);
	std::getline(iss1, line);
	count = stoi(line, nullptr, 10);
	while (std::getline(iss1, line)){
		os << line << endl;
	}
	plane = drawPlaneXY(x, y, -z/2, nDiv, -1);
	std::istringstream iss2 (plane);
	std::getline(iss2, line);
	count += stoi(line, nullptr, 10);
	while (std::getline(iss2, line)){
		os << line << endl;
	}
	plane = drawPlaneXZ(x, y/2, z, nDiv, 1);
	std::istringstream iss3 (plane);
	std::getline(iss3, line);
	count += stoi(line, nullptr, 10);
	while (std::getline(iss3, line)){
		os << line << endl;
	}
	plane = drawPlaneXZ(x, -y/2, z, nDiv, -1);
	std::istringstream iss4 (plane);
	std::getline(iss4, line);
	count += stoi(line, nullptr, 10);
	while (std::getline(iss4, line)){
		os << line << endl;
	}
	plane = drawPlaneYZ(x/2, y, z, nDiv, 1);
	std::istringstream iss5 (plane);
	std::getline(iss5, line);
	count += stoi(line, nullptr, 10);
	while (std::getline(iss5, line)){
		os << line << endl;
	}
	plane = drawPlaneYZ(-x/2, y, z, nDiv, -1);
	std::istringstream iss6 (plane);
	std::getline(iss6, line);
	count += stoi(line, nullptr, 10);
	while (std::getline(iss6, line)){
		os << line << endl;
	}

	return std::to_string(count) + "\n" + os.str();
}
