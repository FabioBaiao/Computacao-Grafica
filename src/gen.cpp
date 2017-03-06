#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <math.h>

using namespace std;

#define UPWARDS 1
#define DOWNWARDS (-1)

string planeXZ (float x, float y, float z, int nDiv, int orient){
	ostringstream os;
	int count = 0;
	float deltaX = x/nDiv;
	float deltaZ = z/nDiv;

	x = x/2;
	z = z/2;

	float xi = x;
	float zi = z;

	for (int i = 1; i <= nDiv; i++){
		for (int j = 1; j <= nDiv; j++){
			if (orient == UPWARDS){
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
	return to_string(count) + "\n" + os.str();
}

string planeXY (float x, float y, float z, int nDiv, int orient){
	ostringstream os;
	int count = 0;
	float deltaX = x/nDiv;
	float deltaY = y/nDiv;

	x = x/2;
	y = y/2;

	float xi = x;
	float yi = -y;

	for (int i = 1; i <= nDiv; i++){
		for (int j = 1; j <= nDiv; j++){
			if (orient == UPWARDS){
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

	return to_string(count) + "\n" + os.str();
}

string planeYZ (float x, float y, float z, int nDiv, int orient){
	ostringstream os;
	int count = 0;
	float deltaY = y/nDiv;
	float deltaZ = z/nDiv;

	y = y/2;
	z = z/2;

	float yi = -y;
	float zi = -z;
	for (int i = 1; i <= nDiv; i++){
		for (int j = 1; j <= nDiv; j++){
			if (orient == UPWARDS){
				os << x << " " << yi << " " << zi << endl;
				os << x << " " << (yi+deltaY) << " " << zi << endl;
				os << x << " " << (yi+deltaY) << " " << (zi+deltaZ) << endl;
				os << x << " " << yi << " " << zi << endl;
				os << x << " " << (yi+deltaY) << " " << (zi+deltaZ) << endl;
				os << x << " " << yi << " " << (zi+deltaZ) << endl;
			} else {
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
	return to_string(count) + "\n" + os.str();
}

string box(float x, float y, float z, int nDiv){
	string plane, line;
	ostringstream os;
	istringstream iss;
	int vertices, count;

	plane = planeXY(x, y, z/2, nDiv, UPWARDS);
	istringstream iss1 (plane);
	getline(iss1, line);
	count = stoi(line, nullptr, 10);
	while (getline(iss1, line)){
		os << line << endl;
	}
	plane = planeXY(x, y, -z/2, nDiv, DOWNWARDS);
	istringstream iss2 (plane);
	getline(iss2, line);
	count += stoi(line, nullptr, 10);
	while (getline(iss2, line)){
		os << line << endl;
	}
	plane = planeXZ(x, y/2, z, nDiv, UPWARDS);
	istringstream iss3 (plane);
	getline(iss3, line);
	count += stoi(line, nullptr, 10);
	while (getline(iss3, line)){
		os << line << endl;
	}
	plane = planeXZ(x, -y/2, z, nDiv, DOWNWARDS);
	istringstream iss4 (plane);
	getline(iss4, line);
	count += stoi(line, nullptr, 10);
	while (getline(iss4, line)){
		os << line << endl;
	}
	plane = planeYZ(x/2, y, z, nDiv, UPWARDS);
	istringstream iss5 (plane);
	getline(iss5, line);
	count += stoi(line, nullptr, 10);
	while (getline(iss5, line)){
		os << line << endl;
	}
	plane = planeYZ(-x/2, y, z, nDiv, DOWNWARDS);
	istringstream iss6 (plane);
	getline(iss6, line);
	count += stoi(line, nullptr, 10);
	while (getline(iss6, line)){
		os << line << endl;
	}

	return to_string(count) + "\n" + os.str();
}

string annulus(float dist, float smj, float smn, int slices) {
	int nPoints = 0;
	float y = 0.0f;
	float alpha, nextAlpha, deltaAlpha;
	ostringstream os;

	deltaAlpha = (2.0f * M_PI) / slices;
	for(int i = 0; i < slices; ++i) {
		float alpha = i * deltaAlpha;
		float nextAlpha = alpha + deltaAlpha;

		os << smj * cosf(alpha) << " " << y << " " << smn * sinf(alpha) << endl;	
		os << (smj+dist) * cosf(nextAlpha) << " " << y << " " << (smn+dist) * sinf(nextAlpha) << endl;	
		os << (smj+dist) * cosf(alpha) << " " << y << " " << (smn+dist) * sinf(alpha) << endl;	
		os << smj * cosf(alpha) << " " << y << " " << smn * sinf(alpha) << endl;	
		os << smj * cosf(nextAlpha) << " " << y << " " << smn * sinf(nextAlpha) << endl;	
		os << (smj+dist) * cosf(nextAlpha) << " " << y << " " << (smn+dist) * sinf(nextAlpha) << endl;	
		nPoints += 6;
	}
	return (to_string(nPoints) + "\n" + os.str());
}

string ellipsoid(float a, float b, float c, float stacks, int slices) {
	int nPoints = 0;
	float deltaAlpha = 2.0f * M_PI / slices;
	float deltaBeta = M_PI / stacks;
	ostringstream os;
		
	for(int i = 0; i < stacks; i++) {
		float beta = i * deltaBeta;
		float nextBeta = beta + deltaBeta;

		for(int j = 0; j < slices; j++) {
			float alpha = j * deltaAlpha;
			float nextAlpha = alpha + deltaAlpha;

			if(i<stacks-1){
				os << a * sinf(beta) * sinf(alpha) << " " << b * cosf(beta) << " " << c * sinf(beta) * cosf(alpha) << endl;
				os << a * sinf(nextBeta) * sinf(alpha) << " " << b * cosf(nextBeta) << " " << c * sinf(nextBeta) * cosf(alpha) << endl;
				os << a * sinf(nextBeta) * sinf(nextAlpha) << " " << b * cosf(nextBeta) << " " << c * sinf(nextBeta) * cosf(nextAlpha) << endl;
				nPoints +=3;
			}
			if(i>0){
				os << a * sinf(beta) * sinf(alpha) << " " << b * cosf(beta) << " " << c * sinf(beta) * cosf(alpha) << endl;
				os << a * sinf(nextBeta) * sinf(nextAlpha) << " " << b * cosf(nextBeta) << " "<< c * sinf(nextBeta) * cosf(nextAlpha) << endl;
				os << a * sinf(beta) * sinf(nextAlpha) << " " << b * cosf(beta) << " " << c * sinf(beta) * cosf(nextAlpha) << endl;
				nPoints +=3;
			}
		}

	}
	return (to_string(nPoints) + "\n" + os.str());
}

string frustum(float baseRadius, float topRadius, float height, int slices, int stacks) {
	int i, j, nPoints = 0;
	float *cosCache, *sinCache;
	float dr, dh, alpha, dAlpha;
	string ret = "";
	ostringstream os;

	cosCache = new float[slices+1];
	sinCache = new float[slices+1];
	alpha = 0.0f;
	dAlpha = 2.0f * M_PI / slices;
	for(i = 0; i < slices; i++) {
		alpha = i * dAlpha;
		cosCache[i] = cosf(alpha);
		sinCache[i] = sinf(alpha);
	}
	cosCache[slices] = cosCache[0]; // cos(2 * M_PI) = cos(0)
	sinCache[slices] = sinCache[0]; // sin(2 * M_PI) = sin(0)

	if(baseRadius > 0.0f) {
		for(i = 0; i < slices; i++) { // lower base
			os << (baseRadius * sinCache[i]) << " 0.0 " << (baseRadius * cosCache[i]) << '\n';
			os << "0.0 0.0 0.0\n";
			os << (baseRadius * sinCache[i+1]) << " 0.0 " << (baseRadius * cosCache[i+1]) << '\n';
			nPoints += 3;
		}
	}
	if(topRadius > 0.0f) {
		for(i = 0; i < slices; ++i) { // top base
			os << (topRadius * sinCache[i]) << ' ' << height << ' ' << (topRadius * cosCache[i]) << '\n';
			os << (topRadius * sinCache[i+1]) << ' ' << height << ' ' << (topRadius * cosCache[i+1]) << '\n';
			os << "0.0 " << height << " 0.0\n";
			nPoints += 3;
		}
	}
	dr = (baseRadius - topRadius) / stacks;
	dh = height / stacks;
	// side of the frustum
	for(i = 0; i < stacks; ++i) {
		float rLow, rHigh;
		float yLow, yHigh;

		rLow = baseRadius - i * dr;
		rHigh = rLow - dr;
		yLow = i * dh;
		yHigh = yLow + dh;
		if(rHigh > 0.0f) {
			for(j = 0; j < slices; ++j) {
				os << (rLow * sinCache[j]) << ' ' << yLow << ' ' << (rLow * cosCache[j]) << '\n';
				os << (rHigh * sinCache[j+1]) << ' ' << yHigh << ' ' << (rHigh * cosCache[j+1]) << '\n';
				os << (rHigh * sinCache[j]) << ' ' << yHigh << ' ' << (rHigh * cosCache[j]) << '\n';
				nPoints += 3;
			}
		}
		if(rLow > 0.0f) {
			for(j = 0; j < slices; ++j) {
				os << (rLow * sinCache[j]) << ' ' << yLow << ' ' << (rLow * cosCache[j]) << '\n';
				os << (rLow * sinCache[j+1]) << ' ' << yLow << ' ' << (rLow * cosCache[j+1]) << '\n';
				os << (rHigh * sinCache[j+1]) << ' ' << yHigh << ' ' << (rHigh * cosCache[j+1]) << '\n';
				nPoints += 3;
			}
		}
	}
	delete[] cosCache; delete[] sinCache;

	return (to_string(nPoints) + "\n" + os.str());
}

string cone(float radius, float height, int slices, int stacks) {
	return frustum(radius, 0.0f, height, slices, stacks);
}

string cylinder(float radius, float height, int slices, int stacks) {
	return frustum(radius, radius, height, slices, stacks);
}

string torus(float innerRadius, float outerRadius, int nsides, int nrings) {
	int i, j, nPoints = 0;
	float deltaPhi, deltaTheta;
	float phi, theta, nextPhi, nextTheta;
	float cosPhi, cosTheta, sinPhi, sinTheta;
	float cosNextPhi, cosNextTheta, sinNextPhi, sinNextTheta;
	string res = "";
	ostringstream os;

	deltaPhi = 2.0f * M_PI / nrings;
	deltaTheta = 2.0f * M_PI / nsides;
	for(i = 0; i < nrings; i++) {
		phi = i * deltaPhi;
		nextPhi = (i + 1) * deltaPhi;
		cosPhi = cosf(phi);
		sinPhi = sinf(phi);
		cosNextPhi = cosf(nextPhi);
		sinNextPhi = sinf(nextPhi);

		for(j = 0; j < nsides; j++, nPoints += 6) {
			float dXZ, nextDXZ;

			theta = j * deltaTheta;
			nextTheta = (j + 1) * deltaTheta;
			cosTheta = cosf(theta);
			sinTheta = sinf(theta);
			cosNextTheta = cosf(nextTheta);
			sinNextTheta = sinf(nextTheta);
			dXZ = outerRadius + innerRadius * cosTheta;
			nextDXZ = outerRadius + innerRadius * cosNextTheta;
			
			os << (dXZ * sinPhi) << ' ' << (innerRadius * sinTheta) << ' ' << (dXZ * cosPhi) << '\n';
			os << (nextDXZ * sinNextPhi) << ' ' << (innerRadius * sinNextTheta) << ' ' << (nextDXZ * cosNextPhi) << '\n';
			os << (nextDXZ * sinPhi) << ' ' << (innerRadius * sinNextTheta) << ' ' << (nextDXZ * cosPhi) << '\n';
			
			os << (dXZ * sinPhi) << ' ' << (innerRadius * sinTheta) << ' ' << (dXZ * cosPhi) << '\n';
			os << (dXZ * sinNextPhi) << ' ' << (innerRadius * sinTheta) << ' ' << (dXZ * cosNextPhi) << '\n';
			os << (nextDXZ * sinNextPhi) << ' ' << (innerRadius * sinNextTheta) << ' ' << (nextDXZ * cosNextPhi) << '\n';
		}
	}
	return (to_string(nPoints) + "\n" + os.str());
}

int annulusGenerator(int argc, char* argv[]){
	float dist, smj, smn;
	int slices;
	
	ofstream outfile;
	dist = atof(argv[0]);
	smj = atof(argv[1]);
	smn = atof(argv[2]);
	slices = atoi(argv[3]);
	if(dist <= 0.0f || smj <= 0.0f || smn <= 0.0f || slices <= 0){
		fputs("Error: All parameters of the annulus must be positive numbers\n", stderr);
		return 1;
	}

	outfile.open(argv[4]);
	if(!outfile.is_open()) {
		perror("ofstream.open");
		return 1;
	}
	outfile << annulus(dist, smj, smn, slices);
	outfile.close();
	return 0;
}

int boxGenerator(int argc, char *argv[]) {
	float xDim, yDim, zDim;
	int divisions;
	ofstream outfile;

	xDim = atof(argv[0]);
	yDim = atof(argv[1]);
	zDim = atof(argv[2]);
	divisions = (argc == 5 ? atoi(argv[3]) : 1);

	if(xDim <= 0.0f || yDim <= 0.0f || zDim <= 0.0f || divisions <= 0) {
		fputs("Error: All parameters of the box must be positive numbers\n", stderr);
		return 1;
	}
	
	outfile.open((argc == 5) ? argv[4] : argv[3]);
	if(!outfile.is_open()) {
		perror("ofstream.open");
		return 1;
	}
	outfile << box(xDim, yDim, zDim, divisions);
	outfile.close();
	return 0;
}

int coneGenerator(int argc, char *argv[]) {
	float radius, height;
	int slices, stacks;
	ofstream outfile;

	radius = atof(argv[0]);
	height = atof(argv[1]);
	slices = atoi(argv[2]);
	stacks = atoi(argv[3]);
	if(radius <= 0.0f || height <= 0.0f || slices <= 0 || stacks <= 0) {
		fputs("Error: All parameters of the cone must be positive numbers\n", stderr);
		return 1;
	}
	outfile.open(argv[4]);
	if(!outfile.is_open()) {
		perror("ofstream.open");
		return 1;
	}
	outfile << cone(radius, height, slices, stacks);
	outfile.close();
	return 0;
}

// Cylinder and cone generators are very similar. Consider reusing code!
int cylinderGenerator(int argc, char *argv[]) {
	float radius, height;
	int slices, stacks;
	ofstream outfile;

	radius = atof(argv[0]);
	height = atof(argv[1]);
	slices = atoi(argv[2]);
	stacks = atoi(argv[3]);
	if(radius <= 0.0f || height <= 0.0f || slices <= 0 || stacks <= 0) {
		fputs("Error: All parameters of the cylinder must be positive numbers\n", stderr);
		return 1;
	}
	outfile.open(argv[4]);
	if(!outfile.is_open()) {
		perror("ofstream.open");
		return 1;
	}
	outfile << cylinder(radius, height, slices, stacks);
	outfile.close();
	return 0;
}

int ellipsoidGenerator(int argc, char *argv[]) {
	float xRadius, yRadius, zRadius;
	int slices, stacks;
	ofstream outfile;

	xRadius = atof(argv[0]);
	yRadius = atof(argv[1]);
	zRadius = atof(argv[2]);
	slices = atoi(argv[3]);
	stacks = atoi(argv[4]);
	if(xRadius <= 0.0f || yRadius <= 0.0f || zRadius <= 0.0f || slices <= 0 || stacks <= 0) {
		fputs("Error: All parameters of the cylinder must be positive numbers\n", stderr); // consider creating a macro for this kind of error message
		return 1;
	}
	outfile.open(argv[5]);
	if(!outfile.is_open()) {
		perror("ofstream.open");
		return 1;
	}
	outfile << ellipsoid(xRadius, yRadius, zRadius, slices, stacks);
	outfile.close();
	return 0;
}

int frustumGenerator(int argc, char *argv[]) {
	float baseRadius, topRadius, height;
	int slices, stacks;
	char *endptr1, *endptr2;
	ofstream outfile;

	baseRadius = strtof(argv[0], &endptr1);
	topRadius = strtof(argv[1], &endptr2);
	height = atof(argv[2]);
	slices = atoi(argv[3]);
	stacks = atoi(argv[4]);
	if(*endptr1 != '\0' || *endptr2 != '\0' || baseRadius < 0.0f || topRadius < 0.0f) {
		fputs("Error: The base radius and the top radius of the frustum must be non-negative numbers\n", stderr);
		return 1;
	}
	if(height <= 0 || slices <= 0 || stacks <= 0) {
		fputs("Error: Height, slices and stacks must all be positive\n", stderr);
		return 1;
	}
	
	outfile.open(argv[5]);
	if(!outfile.is_open()) {
		perror("ofstream.open");
		return 1;
	}
	outfile << frustum(baseRadius, topRadius, height, slices, stacks);
	outfile.close();
	return 0;
}

int planeGenerator(int argc, char *argv[]) {
	float xDim, zDim;
	int divisions;
	ofstream outfile;

	xDim = atof(argv[0]);
	zDim = atof(argv[1]);
	divisions = (argc == 4 ? atoi(argv[2]) : 1);
	if(xDim <= 0.0f || zDim <= 0.0f || divisions <= 0) {
		fputs("Error: All parameters of the plane must be positive numbers\n", stderr);
		return 1;
	}
	(argc == 4) ? outfile.open(argv[3]) : outfile.open(argv[2]);

	if(!outfile.is_open()) {
		perror("ofstream.open");
		return 1;
	}
	outfile << planeXZ(xDim, 0, zDim, divisions, UPWARDS);
	outfile.close();
	return 0;
}

int sphereGenerator(int argc, char *argv[]) {
	float radius;
	int slices, stacks;
	ofstream outfile;

	radius = atof(argv[0]);
	slices = atoi(argv[1]);
	stacks = atoi(argv[2]);
	if(radius <= 0.0f || slices <= 0 || stacks <= 0) {
		fputs("Error: All parameters of the sphere must be positive numbers\n", stderr);
		return 1;
	}

	outfile.open(argv[3]);
	if(!outfile.is_open()) {
		perror("ofstream.open");
		return 1;
	}
	outfile << ellipsoid(radius, radius, radius, slices, stacks);
	outfile.close();
	return 0;
}

int torusGenerator(int argc, char *argv[]) {
	float innerRadius, outerRadius;
	int sides, rings;
	ofstream outfile;

	innerRadius = atof(argv[0]);
	outerRadius = atof(argv[1]);
	sides = atoi(argv[2]);
	rings = atoi(argv[3]);
	if(innerRadius <= 0.0f || outerRadius <= 0.0f || sides <= 0 || rings <= 0) {
		fputs("Error: All parameters of the torus must be positive numbers\n", stderr);
		return 1;
	}
	outfile.open(argv[4]);
	if(!outfile.is_open()) {
		perror("ofstream.open");
		return 1;
	}
	outfile << torus(innerRadius, outerRadius, sides, rings);
	outfile.close();
	return 0;
}



void usage(const char *programName, FILE *stream) {
	fprintf(stream, "Usage: %s primitive parameters outfile\n\n", programName);
	fputs("+-------------+-------------------------------------------+\n"
		  "| primitive   | parameters                                |\n"
		  "+-------------+-------------------------------------------+\n"
		  "| annulus     | distance semimajor semiminor slices       |\n"
		  "| box         | xDim yDim zDim [divisions]                |\n"
		  "| cone        | radius height slices stacks               |\n"
		  "| cylinder    | radius height slices stacks               |\n"
		  "| ellipsoid   | xRadius yRadius zRadius slices stacks     |\n"
		  "| frustum     | baseRadius topRadius height slices stacks |\n"
		  "| plane       | xDim zDim [divisions]                     |\n"
		  "| sphere      | radius slices stacks                      |\n"
		  "| torus       | innerRadius outerRadius sides rings       |\n"
		  "+-------------+-------------------------------------------+\n"
		  , stream
	);
}

int main(int argc, char *argv[]) {
	int rval;

	if(argc == 1){
		usage(argv[0], stderr);
		return 1;
	} 

	string primitive(argv[1]); 
	if(primitive == "annulus" && argc == 7) {
		rval = annulusGenerator(argc - 2, &argv[2]);
	}
	else if(primitive == "box" && (argc == 6 || argc == 7)) {
		rval = boxGenerator(argc - 2, &argv[2]);
	}
	else if(primitive == "cone" && argc == 7){
		rval = coneGenerator(argc - 2, &argv[2]);
	}
	else if(primitive == "cylinder" && argc == 7){
		rval = cylinderGenerator(argc - 2, &argv[2]);
	}
	else if(primitive == "ellipsoid" && argc == 8) {
		rval = ellipsoidGenerator(argc - 2, &argv[2]);
	}
	else if(primitive == "frustum" && argc) {
		rval = frustumGenerator(argc - 2, &argv[2]);
	}
	else if(primitive == "plane" && (argc == 5 || argc == 6)) {
		rval = planeGenerator(argc - 2, &argv[2]);
	}
	else if(primitive == "sphere" && argc == 6){
		rval = sphereGenerator(argc - 2, &argv[2]);
	}
	else if(primitive == "torus" && argc == 7){
		rval = torusGenerator(argc - 2, &argv[2]);
	}
	else if(primitive == "--help") {
		usage(argv[0], stdout);
		rval = 0;
	} else {
		usage(argv[0], stderr);
		rval = 1;
	}
	return rval;
}
