#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <string.h>
#include <stdlib.h>


using namespace std;

#define UPWARDS 1
#define DOWNWARDS (-1)

#define BUFF_SIZE 1024

typedef struct point{
	float x;
	float y;
	float z;
}point;

point* cpoints;		// control points
int** indexes;		// indexes dos pontos dos varios patches
int patches;		// numero de patches
int ncpoints;		// numero de control points

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
				os << xi << " " << y << " " << zi << '\n';
				os << xi << " " << y << " " << (zi-deltaZ) << '\n';
				os << (xi-deltaX) << " " << y << " " << (zi-deltaZ) << '\n';

				os << xi << " " << y << " " << zi << '\n';
				os << (xi-deltaX) << " " << y << " " << (zi-deltaZ) << '\n';
				os << (xi-deltaX) << " " << y << " " << zi << '\n';
			}
			else {
				os << xi << " " << y << " " << zi << '\n';
				os << (xi-deltaX) << " " << y << " " << zi << '\n';
				os << (xi-deltaX) << " " << y << " " << (zi-deltaZ) << '\n';

				os << xi << " " << y << " " << zi << '\n';
				os << (xi-deltaX) << " " << y << " " << (zi-deltaZ) << '\n';
				os << xi << " " << y << " " << (zi-deltaZ) << '\n';
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
				os << xi << " " << yi << " " << z << '\n';
				os << xi << " " << (yi+deltaY) << " " << z << '\n';
				os << (xi-deltaX) << " " << (yi+deltaY) << " " << z << '\n'; 
				os << xi << " " << yi << " " << z << '\n';
				os << (xi-deltaX) << " " << (yi+deltaY) << " " << z << '\n';
				os << (xi-deltaX) << " " << yi << " " << z << '\n';
			}
			else {
				os << xi << " " << yi << " " << z << '\n';
				os << (xi-deltaX) << " " << yi << " " << z << '\n';
				os << (xi-deltaX) << " " << (yi+deltaY) << " " << z << '\n'; 
				os << xi << " " << yi << " " << z << '\n';
				os << (xi-deltaX) << " " << (yi+deltaY) << " " << z << '\n';
				os << xi << " " << (yi+deltaY) << " " << z << '\n';
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
				os << x << " " << yi << " " << zi << '\n';
				os << x << " " << (yi+deltaY) << " " << zi << '\n';
				os << x << " " << (yi+deltaY) << " " << (zi+deltaZ) << '\n';
				os << x << " " << yi << " " << zi << '\n';
				os << x << " " << (yi+deltaY) << " " << (zi+deltaZ) << '\n';
				os << x << " " << yi << " " << (zi+deltaZ) << '\n';
			} else {
				os << x << " " << yi << " " << zi << '\n';
				os << x << " " << yi << " " << (zi+deltaZ) << '\n';
				os << x << " " << (yi+deltaY) << " " << (zi+deltaZ) << '\n';
				os << x << " " << yi << " " << zi << '\n';
				os << x << " " << (yi+deltaY) << " " << (zi+deltaZ) << '\n';
				os << x << " " << (yi+deltaY) << " " << zi << '\n';
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
	int count;

	plane = planeXY(x, y, z/2, nDiv, UPWARDS);
	istringstream iss1 (plane);
	getline(iss1, line);
	count = stoi(line, nullptr, 10);
	while (getline(iss1, line)){
		os << line << '\n';
	}
	plane = planeXY(x, y, -z/2, nDiv, DOWNWARDS);
	istringstream iss2 (plane);
	getline(iss2, line);
	count += stoi(line, nullptr, 10);
	while (getline(iss2, line)){
		os << line << '\n';
	}
	plane = planeXZ(x, y/2, z, nDiv, UPWARDS);
	istringstream iss3 (plane);
	getline(iss3, line);
	count += stoi(line, nullptr, 10);
	while (getline(iss3, line)){
		os << line << '\n';
	}
	plane = planeXZ(x, -y/2, z, nDiv, DOWNWARDS);
	istringstream iss4 (plane);
	getline(iss4, line);
	count += stoi(line, nullptr, 10);
	while (getline(iss4, line)){
		os << line << '\n';
	}
	plane = planeYZ(x/2, y, z, nDiv, UPWARDS);
	istringstream iss5 (plane);
	getline(iss5, line);
	count += stoi(line, nullptr, 10);
	while (getline(iss5, line)){
		os << line << '\n';
	}
	plane = planeYZ(-x/2, y, z, nDiv, DOWNWARDS);
	istringstream iss6 (plane);
	getline(iss6, line);
	count += stoi(line, nullptr, 10);
	while (getline(iss6, line)){
		os << line << '\n';
	}

	return to_string(count) + "\n" + os.str();
}

string annulus(float dist, float smj, float smn, int slices) {
	int nPoints = 0;
	float deltaAlpha;
	ostringstream os;

	deltaAlpha = (2.0f * M_PI) / slices;
	for(int i = 0; i < slices; ++i) {
		float alpha = i * deltaAlpha;
		float nextAlpha = alpha + deltaAlpha;

		os << smj * cosf(alpha) << " 0.0 " << smn * sinf(alpha) << '\n';
		os << (smj+dist) * cosf(nextAlpha) << " 0.0 " << (smn+dist) * sinf(nextAlpha) << '\n';
		os << (smj+dist) * cosf(alpha) << " 0.0 " << (smn+dist) * sinf(alpha) << '\n';
		os << smj * cosf(alpha) << " 0.0 " << smn * sinf(alpha) << '\n';
		os << smj * cosf(nextAlpha) << " 0.0 " << smn * sinf(nextAlpha) << '\n';
		os << (smj+dist) * cosf(nextAlpha) << " 0.0 " << (smn+dist) * sinf(nextAlpha) << '\n';
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
				os << a * sinf(beta) * sinf(alpha) << " " << b * cosf(beta) << " " << c * sinf(beta) * cosf(alpha) << '\n';
				os << a * sinf(nextBeta) * sinf(alpha) << " " << b * cosf(nextBeta) << " " << c * sinf(nextBeta) * cosf(alpha) << '\n';
				os << a * sinf(nextBeta) * sinf(nextAlpha) << " " << b * cosf(nextBeta) << " " << c * sinf(nextBeta) * cosf(nextAlpha) << '\n';
				nPoints += 3;
			}
			if(i>0){
				os << a * sinf(beta) * sinf(alpha) << " " << b * cosf(beta) << " " << c * sinf(beta) * cosf(alpha) << '\n';
				os << a * sinf(nextBeta) * sinf(nextAlpha) << " " << b * cosf(nextBeta) << " "<< c * sinf(nextBeta) * cosf(nextAlpha) << '\n';
				os << a * sinf(beta) * sinf(nextAlpha) << " " << b * cosf(beta) << " " << c * sinf(beta) * cosf(nextAlpha) << '\n';
				nPoints += 3;
			}
		}

	}
	return (to_string(nPoints) + "\n" + os.str());
}

string frustum(float baseRadius, float topRadius, float height, int slices, int stacks) {
	int i, j, nPoints = 0;
	float *cosCache, *sinCache;
	float dr, dh, alpha, dAlpha;
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

void cross(float *a, float *b, float *res) {

	res[0] = a[1]*b[2] - a[2]*b[1];
	res[1] = a[2]*b[0] - a[0]*b[2];
	res[2] = a[0]*b[1] - a[1]*b[0];
}

float length(float *v) {

	float res = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
	return res;

}

void normalize(float *a) {

	float l = length(a);
	a[0] = a[0]/l;
	a[1] = a[1]/l;
	a[2] = a[2]/l;
}

void multMatrixVector(float *m, float *v, float *res) {

	for (int j = 0; j < 4; ++j) {
		res[j] = 0;
		for (int k = 0; k < 4; ++k) {
			res[j] += v[k] * m[j * 4 + k];
		}
	}
}

void multVectorMatrix(float *v, float *m, float *res) {
	for(int i = 0; i < 4; ++i){
		res[i] = 0;
		for(int j = 0; j < 4; ++j){
			res[i] += v[j] * m[j*4+i]; 
		}
	}
}

void multMatrixMatrix(float *m1, float *m2, float *res){
	for(int i = 0; i < 4; ++i){
		for(int j = 0; j < 4; ++j){
			res[i*4+j] = 0.0f;
			for(int k = 0; k < 4; ++k)
				res[i*4+j] += m1[i*4+k] * m2[k*4+j];
		}
	}
}

void getBezierPatchPoint(float u, float v, point* pv, float *res, float* p_normal){
	float d_u[3];
	float d_v[3];
	/* Setup */	
	float m[4][4] = {	{-1.0f, 3.0f, -3.0f, 1.0f},
						{ 3.0f,-6.0f,  3.0f, 0.0f},
						{-3.0f, 3.0f,  0.0f, 0.0f},
						{ 1.0f, 0.0f,  0.0f, 0.0f}};

	/* a transposta de m é igual a m */

	float Px[4][4] = { 	{ pv[0].x, pv[1].x, pv[2].x, pv[3].x},
						{ pv[4].x, pv[5].x, pv[6].x, pv[7].x},
						{ pv[8].x, pv[9].x, pv[10].x, pv[11].x},
						{ pv[12].x, pv[13].x, pv[14].x, pv[15].x}};
	
	float Py[4][4] = { 	{ pv[0].y, pv[1].y, pv[2].y, pv[3].y},
						{ pv[4].y, pv[5].y, pv[6].y, pv[7].y},
						{ pv[8].y, pv[9].y, pv[10].y, pv[11].y},
						{ pv[12].y, pv[13].y, pv[14].y, pv[15].y}};
	
	float Pz[4][4] = { 	{ pv[0].z, pv[1].z, pv[2].z, pv[3].z},
						{ pv[4].z, pv[5].z, pv[6].z, pv[7].z},
						{ pv[8].z, pv[9].z, pv[10].z, pv[11].z},
						{ pv[12].z, pv[13].z, pv[14].z, pv[15].z}};
	
	float U[4] = {u*u*u, u*u, u, 1};
	float UD[4]= {3*u*u, 2*u, 1, 0};
	float V[4] = {v*v*v, v*v, v, 1};
	float VD[4]= {3*v*v, 2*v, 1, 0};

	// solution 1: doesnt work and i dont know why
	float MdV[4];
	float MV[4];
	multMatrixVector((float*)m, V, MV);
	multMatrixVector((float*)m, VD, MdV);

	float dUM[4];
	float UM[4];
	multVectorMatrix(U, (float*)m, UM);
	multVectorMatrix(UD, (float*)m, dUM);
	
	float UMP[3][4];
	multVectorMatrix(UM, (float*) Px, UMP[0]);	
	multVectorMatrix(UM, (float*) Py, UMP[1]);	
	multVectorMatrix(UM, (float*) Pz, UMP[2]);	

	float dUMP[3][4];
	multVectorMatrix(dUM, (float*) Px, dUMP[0]);	
	multVectorMatrix(dUM, (float*) Py, dUMP[1]);	
	multVectorMatrix(dUM, (float*) Pz, dUMP[2]);	
	
	for(int j = 0; j < 3; j++)
	{
		res[j] = 0.0f;
		d_u[j] = 0.0f;
		d_v[j] = 0.0f;
		for(int i = 0; i < 4; i++ )
		{
			res[j] += MV[i] * UMP[j][i]; 
			d_u[j] += MV[i] * dUMP[j][i];
			d_v[j] += MdV[i] * UMP[j][i];
		}
	}
	
	normalize(d_u);
	normalize(d_v);
	cross(d_v, d_u, p_normal);
}


void bezierPatchGenerator(char* outfile)
{
	point pv[16];
	int divs = 16; // change this to change the tesselation level
	int p = 0;
	ofstream out;
	ostringstream os;
	
	for(int i = 0; i < patches; i++)
	{
		for(int j = 0; j < 16; j++)
			pv[j] = cpoints[indexes[i][j]];
		
		for(int u = 0; u < divs; u++)
		{
			float resP1[3], p_normal1[3];
			float resP2[3], p_normal2[3];
			float resP3[3], p_normal3[3];
			float resP4[3], p_normal4[3];
			
			for(int v = 0; v < divs; v++)
			{
				
				getBezierPatchPoint(u/(float)divs, v/(float)divs, pv, resP1,p_normal1 );
				getBezierPatchPoint((u+1)/(float)divs, v/(float)divs, pv, resP2, p_normal2 );
				getBezierPatchPoint(u/(float)divs, (v+1)/(float)divs, pv, resP3, p_normal3);
				getBezierPatchPoint((u+1)/(float)divs, (v+1)/(float)divs, pv, resP4,p_normal4);
				
				os << resP1[0] << ' ' << resP1[1] << ' ' << resP1[2] << '\n';
				os << resP3[0] << ' ' << resP3[1] << ' ' << resP3[2] << '\n';
				os << resP4[0] << ' ' << resP4[1] << ' ' << resP4[2] << '\n';
				
				os << resP2[0] << ' ' << resP2[1] << ' ' << resP2[2] << '\n';
				os << resP1[0] << ' ' << resP1[1] << ' ' << resP1[2] << '\n';
				os << resP4[0] << ' ' << resP4[1] << ' ' << resP4[2] << '\n';
				
				p += 6;
			}
		}
	}
	out.open(outfile);
	out << to_string(p) + "\n" + os.str();
	out.close();	
}

int bezierPatchParser(char* patch)
{
	int l_index;
	int i, j;
	char line[BUFF_SIZE];
	FILE* f = fopen(patch, "r");
	if(!f)
		return -1;

	fscanf(f, "%d\n", &l_index);
	patches = l_index;
	printf("%d\n", patches);
	indexes = (int**) malloc(sizeof(int*)*l_index);
	if(!indexes)
		return -1;

	for(i = 0; i < l_index; i++)
	{
		indexes[i] = (int*) malloc(sizeof(int) * 16);
		if(!indexes[i])
			return -1;
		
		memset(line, 0, BUFF_SIZE);
		fgets(line, BUFF_SIZE, f);
		char* ind = NULL;
		for(j = 0, ind = strtok(line,", "); ind && j < 16; ind = strtok(NULL, ", "), j++)
			indexes[i][j] = atoi(ind);
	}
	
	fscanf(f, "%d\n", &ncpoints);

	cpoints = (point*) malloc(sizeof(point) * ncpoints);
	if(!cpoints)
		return -1;
	
	for(i = 0; i < ncpoints; i++)
	{
		memset(line, 0, BUFF_SIZE);
		fgets(line,BUFF_SIZE, f);
		cpoints[i].x = atof(strtok(line, ", "));
		cpoints[i].y = atof(strtok(NULL, ", "));
		cpoints[i].z = atof(strtok(NULL, ", "));
	}

	fclose(f);
	return 0;
}

void usage(const char *programName, FILE *stream) {
	fprintf(stream, "Usage: %s primitive parameters outfile\n\n", programName);
	fputs("+-------------+-------------------------------------------+\n"
		  "| primitive   | parameters                                |\n"
		  "+-------------+-------------------------------------------+\n"
		  "| annulus     | thickness semimajor semiminor slices      |\n"
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
	else if(primitive == "frustum" && argc == 8) {
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
	else if(primitive == "patch" && argc == 4)
	{
		rval = bezierPatchParser(argv[2]);
		if(!rval)
			bezierPatchGenerator(argv[3]);
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
