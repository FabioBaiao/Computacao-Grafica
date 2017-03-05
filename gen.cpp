#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <math.h>

using namespace std;

std::string sphere(float radius, float slices, int stacks) {

	std:string ret = "";
	int count = 0;
	int half_stacks = (int) stacks / 2;
	float b_step = M_PI / stacks;
	float a_step = 2*M_PI / slices;
		
	for(float i = M_PI/2; i > -M_PI/2; i -= b_step)
	{
//		float beta = i * (M_PI) / stacks;
		float delta_b = i - (M_PI) / stacks;
		
		for(int j = 0; j < slices; j++)
		{
			float alfa = j * (2*M_PI)/slices;
			float delta_a = alfa + (2*M_PI)/slices;

			ret += std::to_string(radius*cos(i)*sin(alfa)) + " " + 
				   std::to_string(radius*sin(i)) + " " + 
				   std::to_string(radius*cos(i)*cos(alfa)) + "\n";
//			glVertex3f(radius*cos(i)*sin(alfa), radius*sin(i), radius*cos(i)*cos(alfa));
			
			ret += std::to_string(radius*cos(delta_a)*sin(alfa)) + " " + 
				   std::to_string(radius*sin(delta_a)) + " " + 
				   std::to_string(radius*cos(delta_a)*cos(alfa)) + "\n";
//			glVertex3f(radius*cos(delta_b)*sin(alfa), radius*sin(delta_b), radius*cos(delta_b)*cos(alfa));
			
			ret += std::to_string(radius*cos(delta_a)*sin(delta_a)) + " " + 
				   std::to_string(radius*sin(delta_b)) + " " + 
				   std::to_string(radius*cos(delta_b)*cos(delta_a)) + "\n";
//			glVertex3f(radius*cos(delta_b)*sin(delta_a), radius*sin(delta_b), radius*cos(delta_b)*cos(delta_a));
			
			ret += std::to_string(radius*cos(i)*sin(alfa)) + " " + 
				   std::to_string(radius*sin(i)) + " " + 
				   std::to_string(radius*cos(i)*cos(alfa)) + "\n";
//			glVertex3f(radius*cos(i)*sin(alfa), radius*sin(i), radius*cos(i)*cos(alfa));
			
			ret += std::to_string(radius*cos(delta_b)*sin(delta_a)) + " " + 
				   std::to_string(radius*sin(delta_b)) + " " + 
				   std::to_string(radius*cos(delta_b)*cos(delta_a)) + "\n";
//			glVertex3f(radius*cos(delta_b)*sin(delta_a), radius*sin(delta_b), radius*cos(delta_b)*cos(delta_a));
			
			ret += std::to_string(radius*cos(i)*sin(delta_a)) + " " + 
				   std::to_string(radius*sin(i)) + " " + 
				   std::to_string(radius*cos(i)*cos(delta_a)) + "\n";
//			glVertex3f(radius*cos(i)*sin(delta_a), radius*sin(i), radius*cos(i)*cos(delta_a));
	
			count+=6;
		}

	}

	return std::to_string(count) + "\n" + ret; 
}

std::string frustum(float baseRadius, float topRadius, float height, int slices, int stacks) {
	int i, j, nPoints = 0;
	float *cosCache, *sinCache;
	float dr, dh, alpha, dAlpha;
	std::string ret = "";
	std::ostringstream os;

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

	return (std::to_string(nPoints) + os.str());
}

std::string cone(float radius, float height, int slices, int stacks) {
	return frustum(radius, 0.0f, height, slices, stacks);
}

std::string cylinder(float radius, float height, int slices, int stacks) {
	return frustum(radius, radius, height, slices, stacks);
}

std::string torus(float innerRadius, float outerRadius, int nsides, int nrings) {
	int i, j, nPoints = 0;
	float deltaPhi, deltaTheta;
	float phi, theta, nextPhi, nextTheta;
	float cosPhi, cosTheta, sinPhi, sinTheta;
	float cosNextPhi, cosNextTheta, sinNextPhi, sinNextTheta;
	std::string res = "";
	std::ostringstream os;

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
	return (std::to_string(nPoints) + os.str());
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
	// outfile << box(xDim, yDim, zDim, divisions);
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
	// outfile << ellipsoid(xRadius, yRadius, zRadius, slices, stacks);
	outfile.close();
	return 0;
}

int frustumGenerator(int argc, char *argv[]) {
	float baseRadius, topRadius, height;
	int slices, stacks;
	ofstream outfile;

	baseRadius = atof(argv[0]);
	topRadius = atof(argv[1]);
	height = atof(argv[2]);
	slices = atoi(argv[3]);
	stacks = atoi(argv[4]);
	if(baseRadius < 0.0f || topRadius < 0.0f) {
		fputs("Error: The base radius and the top radius must both be non-negative\n", stderr);
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
	float xDim, yDim;
	int divisions;
	ofstream outfile;

	xDim = atof(argv[0]);
	yDim = atof(argv[1]);
	divisions = (argc == 4 ? atoi(argv[2]) : 1);
	if(xDim <= 0.0f || yDim <= 0.0f || divisions <= 0) {
		fputs("Error: All parameters of the plane must be positive numbers\n", stderr);
		return 1;
	}
	(argc == 4) ? outfile.open(argv[3]) : outfile.open(argv[2]);

	if(!outfile.is_open()) {
		perror("ofstream.open");
		return 1;
	}
	// outfile << plane(xDim, yDim, divisions);
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
	outfile << sphere(radius, slices, stacks);
	outfile.close();
	return 0;
}

int torusGenerator(int argc, char *argv[]) {
	float innerRadius, outerRadius;
	int stacks, sides, rings;
	ofstream outfile;

	innerRadius = atof(argv[0]);
	outerRadius = atof(argv[1]);
	sides = atoi(argv[2]);
	stacks = atoi(argv[3]);
	if(innerRadius <= 0.0f || outerRadius <= 0.0f || sides <= 0 || rings <= 0) {
		fputs("Error: All parameters of the torus must be positive numbers\n", stderr);
		return 1;
	}
	outfile.open(argv[3]);
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
		  "| box         | xDim yDim zDim [divisions]                |\n"
		  "| cone        | radius height slices stacks               |\n"
		  "| cylinder    | radius height slices stacks               |\n"
		  "| ellipsoid   | xRadius yRadius zRadius slices stacks     |\n"
		  "| frustum     | baseRadius topRadius height slices stacks |\n"
		  "| plane       | xDim yDim [divisions]                     |\n"
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

	std::string primitive(argv[1]); 
	//std::vector<Point> points;
	if(primitive == "box" && (argc == 6 || argc == 7)) {
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


