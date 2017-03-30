#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <ctype.h>
#define _USE_MATH_DEFINES
#include <math.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "geoTransform.h"
#include "group.h"
#include "point.h"
#include "rotation.h"
#include "scale.h"
#include "tinyxml2.h"
#include "translation.h"

using namespace std;
using namespace tinyxml2;

typedef vector<point> Model;

// Camera control
float r = 10.0f;
float alpha;
float beta;

float pitch = 0.0f, yaw = 0.0f;

float Px = 2.0f, Py = 2.0f, Pz = 2.0f;
float lookX = Px + cos(pitch) * sin(yaw);
float lookY = Py + sin(pitch);
float lookZ = Pz + cos(pitch) * cos(yaw);
// Polygon Mode
GLenum mode;

vector<group> groups;
map<string, Model> models;

// directory of the read file
string directory; 

// TODO!!
// need to correct to accept full paths
string directoryOfFile(const string& fname) {
	size_t pos = fname.find_last_of("\\/");

	return (string::npos == pos)? "" : fname.substr(0, pos+1);
}

float randFloat() {
	 return static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
}

void changeSize(int w, int h) {
	if(h == 0)
		h = 1;
	
	float ratio = w * 1.0 / h;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0, 0, w, h);
	gluPerspective(45.0f ,ratio, 1.0f ,1000.0f);
	glMatrixMode(GL_MODELVIEW);
}

void print_matrix(float m[], int I, int J) {
	for(int i  = 0; i < I; i++){
		for(int j  = 0; j < J; j++){
			cout << m[i + j * I] << " ";
		}
		cout << endl;
	}
}

void drawGroup(group g) {
	glPushMatrix();
	
	for(geoTransform* t : g.transforms) {
		t->apply();
	}
	int i = 0;
	for(string model : g.models) {
		color c = g.modelsColor[i];
		
		glBegin(GL_TRIANGLES);
		glColor3f( c.r, c.g, c.b);
		for(auto p : models[model]) {
			glVertex3f(p.x, p.y, p.z);
		}
		glEnd();
		i++;
	}
	for(auto gr : g.childGroups) {
		drawGroup(gr);
	}
	glPopMatrix();
}

void renderScene(void) {
    GLenum modes[] = {GL_FILL, GL_LINE, GL_POINT};
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPolygonMode(GL_FRONT, modes[mode]);
	glLoadIdentity();
/*	gluLookAt(r*cosf(beta)*cosf(alpha), r*sinf(beta), r*cosf(beta)*sinf(alpha), 
		      0.0,0.0,0.0,
			  0.0f,1.0f,0.0f);
*/
	gluLookAt(Px, Py, Pz, 
	          lookX,lookY,lookZ,
		      0.0f,1.0f,0.0f);
	for(auto g : groups) {
		drawGroup(g);
	}
	glutSwapBuffers();
}

void lookAt(float alpha, float beta){
	lookX = Px + cos(beta) * sin(alpha);
	lookY = Py + sin(beta);
	lookZ = Pz + cos(beta) * cos(alpha);
}

float magnitude(float vx, float vy, float vz){
	return sqrtf((vx*vx) + (vy*vy) + (vz*vz));
}

float normalize(float* vx, float* vy, float* vz){
	float mag = magnitude(*vx, *vy, *vz);
	*vx = (*vx) / mag;
	*vy = (*vy) / mag;
	*vz = (*vz) / mag;
}

void processKeys(unsigned char c, int xx, int yy) {
// put code to process regular keys in here
	
	/* */
	float k = 0.5f;
	
	/*up vector*/
	float upX = 0.0f, upY = 1.0f, upZ = 0.0f;

	/* forward vector */
	float dX = lookX - Px;
	float dY = lookY - Py;
	float dZ = lookZ - Pz;

	normalize(&dX, &dY, &dZ);

	float rX, rY, rZ;

	switch(toupper(c)){
		case 27: exit(0);
		case 'W': 
			Px += k*dX;
			Py += k*dY;
			Pz += k*dZ;

			lookX += k*dX;
			lookY += k*dY;
			lookZ += k*dZ;
			break;
		case 'S': 
			Px -= k*dX;
			Py -= k*dY;
			Pz -= k*dZ;

			lookX -= k*dX;
			lookY -= k*dY;
			lookZ -= k*dZ;
			break;
		case 'A': 
			/* cross product: up x forwardV */
			rX = (upY * dZ) - (upZ * dY);
			rY = (upZ * dX) - (upX * dZ);
			rZ = (upX * dY) - (upY * dX);

			Px += k*rX;
			Py += k*rY;
			Pz += k*rZ;
	
			lookX += k*rX;
			lookY += k*rY;
			lookZ += k*rZ;
			break;
		case 'D': 
			/* cross product: up x forwardV */
			rX = (dY * upZ) - (dZ * upY);
			rY = (dZ * upX) - (dX * upZ);
			rZ = (dX * upY) - (dY * upX);

			Px += k*rX;
			Py += k*rY;
			Pz += k*rZ;
	
			lookX += k*rX;
			lookY += k*rY;
			lookZ += k*rZ;
			break;
		case 'M': // More radius
			r += 0.2f;
			break;
		case 'L': // Less radius
			r -= 0.2f;
			if(r < 0.2f)
				r = 0.2f;
			break;
		case 'C':
			mode = (mode + 1) % 3;
			break;
	}
	glutPostRedisplay();
}

void processSpecialKeys(int key, int xx, int yy) {
// put code to process special keys in here
	switch(key){
		case GLUT_KEY_UP:
			pitch += 0.01f;
			if(pitch > 1.5f)
				pitch = 1.5f;
			break;
		case GLUT_KEY_LEFT:
			yaw += 0.01f;
			break;
		case GLUT_KEY_DOWN:
			pitch -= 0.01f;
			if(pitch < -1.5f)
				pitch = -1.5;
			break;
		case GLUT_KEY_RIGHT:
			yaw -= 0.01f;
			break;
	}
	lookAt(yaw, pitch);
	glutPostRedisplay();
}

group parseGroup(XMLElement *gr) {
	group g;
	XMLElement *child = gr->FirstChildElement();
	for( ; child; child = child->NextSiblingElement()) {
		string type = string(child->Name()); 
		
		if(type == "translate"){
			float x, y, z;
			x = y = z = 0.0f;

			child->QueryFloatAttribute("X", &x);
			child->QueryFloatAttribute("Y", &y);
			child->QueryFloatAttribute("Z", &z);
			
			g.transforms.push_back(new translation(x,y,z));
		} else if(type == "rotate") {
			float angle, axisX, axisY, axisZ;
			angle = axisX = axisY = axisZ = 0.0f;

			child->QueryFloatAttribute("angle", &angle);
			child->QueryFloatAttribute("axisX", &axisX);
			child->QueryFloatAttribute("axisY", &axisY);
			child->QueryFloatAttribute("axisZ", &axisZ);

			g.transforms.push_back(new rotation(angle, axisX, axisY, axisZ));
		} else if(type == "scale") {
			float x, y, z;
			x = y = z = 1.0f;

			child->QueryFloatAttribute("X", &x);
			child->QueryFloatAttribute("Y", &y);
			child->QueryFloatAttribute("Z", &z);

			g.transforms.push_back(new scale(x,y,z));
		} else if(type == "group") {
			group g_child = parseGroup(child);
			g.childGroups.push_back(g_child);
		} else if(type == "models") {
			XMLElement *model = child->FirstChildElement("model");
			
			for(; model; model=model->NextSiblingElement()){
				const char* filename= model->Attribute("file");
				
				if(filename != NULL) {
					string fName = string(filename);
					int n_vertex;

					ifstream file(directory + filename);

					if(!file) {
						cerr << "The file \"" << filename << "\" was not found.\n";
					}

					float rr, gg, bb;
					int r_r, r_g, r_b;
					rr = gg = bb = 0.0;

					r_r = model->QueryFloatAttribute("R", &rr);
					r_g = model->QueryFloatAttribute("G", &gg);
					r_b = model->QueryFloatAttribute("B", &bb);

					if(r_r != XML_SUCCESS && r_g != XML_SUCCESS && r_b != XML_SUCCESS){
						// the color defaults to white if not specified
						rr = gg = bb = 1;
					}
					color c(rr,gg,bb);

					g.models.push_back(fName);
					g.modelsColor.push_back(c);
					if(models.find(fName) != models.end()){
						// if the model file was already read
						continue;
					}

					// reads the number of vertices from the file
					file >> n_vertex; 
					Model model_read;

					for(int i = 0; i < n_vertex; i++){
						float color_r, color_g, color_b;
						float px, py, pz;
						file >> px;
						file >> py;
						file >> pz;
						point p(px, py, pz);
						model_read.push_back(p);
					}
					file.close();
					models[fName] = model_read;
				}
			}
		}
	}
	return g;
}

//We assume that the .xml and .3d files passed are correct.
int main(int argc, char **argv) {
	if(argc != 2) {
		cerr << "Usage: " << argv[0] << " config_file\n";
		return 1;
	}
	XMLDocument doc;
	XMLError loadOkay = doc.LoadFile(argv[1]);

	if(loadOkay != XML_SUCCESS) { 
		perror("LoadFile");
		cerr << "Error loading file '" << argv[1] << "'.\n";
		return 1;
	}

	directory = directoryOfFile(argv[1]);
	XMLElement* scene = doc.FirstChildElement("scene");
	if(scene == NULL) {
		cerr << "Error: The first XML element must be \"scene\"" << endl;
		return 1;
	}

	XMLElement* gr = scene->FirstChildElement("group");
	for(; gr; gr = gr->NextSiblingElement()){
		group g = parseGroup(gr);
		groups.push_back(g);
	}

	// init GLUT and the window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH|GLUT_DOUBLE|GLUT_RGBA);
	glutInitWindowPosition(100,100);
	glutInitWindowSize(800,800);
	glutCreateWindow("Pratical Assignment");

	// Required callback registry 
	glutDisplayFunc(renderScene);
	glutReshapeFunc(changeSize);
	glutKeyboardFunc(processKeys);
	glutSpecialFunc(processSpecialKeys);

	//  OpenGL settings
	glEnable(GL_DEPTH_TEST);
	glMatrixMode(GL_MODELVIEW);
	glEnable(GL_CULL_FACE);

	glutMainLoop();

	return 1;
}
