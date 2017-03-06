#define TIXML_USE_STL 
#include "tinyxml/tinyxml.h"
#include "point.h"
#include "triangle.h"
#include <vector>
#include <fstream>
#include <string>
#define _USE_MATH_DEFINES
#include <math.h>
#include <ctype.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

using namespace std;

// Camera control
float r = 10.0f;
float alpha;
float beta;

// Polygon Mode
GLenum mode;

// a figure is a set of triangles 
typedef vector<triangle> figure;

// Structure to save figures to draw
vector<figure> figures;


string directory_of_file(const string& fname)
{
	size_t pos = fname.find_last_of("\\/");
	return (string::npos == pos)? "" : fname.substr(0, pos+1);
}

float randFloat() {
	 return static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
}

void changeSize(int w, int h) {

	// Prevent a divide by zero, when window is too short
	// (you cant make a window with zero width).
	if(h == 0)
		h = 1;

	// compute window's aspect ratio 
	float ratio = w * 1.0 / h;

	// Set the projection matrix as current
	glMatrixMode(GL_PROJECTION);
	// Load Identity Matrix
	glLoadIdentity();
	
	// Set the viewport to be the entire window
	glViewport(0, 0, w, h);

	// Set perspective
	gluPerspective(45.0f ,ratio, 1.0f ,1000.0f);

	// return to the model view matrix mode
	glMatrixMode(GL_MODELVIEW);
}

void drawTriangle(triangle t){
	glColor3f(t.color_r, t.color_g, t.color_b);
	glVertex3f(t.p1.x, t.p1.y, t.p1.z);
	glVertex3f(t.p2.x, t.p2.y, t.p2.z);
	glVertex3f(t.p3.x, t.p3.y, t.p3.z);
}

void renderScene(void) {
	// clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	// set the camera
	glLoadIdentity();
	gluLookAt(r*cosf(beta)*cosf(alpha), r*sinf(beta), r*cosf(beta)*sinf(alpha), 
		      0.0,0.0,0.0,
			  0.0f,1.0f,0.0f);
	
	for(auto fig:figures) {
		int size = fig.size();
        GLenum modes[] = {GL_FILL, GL_LINE, GL_POINT};
        glPolygonMode(GL_FRONT, modes[mode]);

		glBegin(GL_TRIANGLES);
		for(int i = 0; i < size; i++){
			drawTriangle(fig[i]);
		}
		glEnd();
	}
	// End of frame
	glutSwapBuffers();
}


void processKeys(unsigned char c, int xx, int yy) {
// put code to process regular keys in here
	switch(toupper(c)){
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
			beta += 0.1f;
			if(beta > 1.5f)
				beta = 1.5f;
			break;
		case GLUT_KEY_LEFT:
			alpha += 0.1f;
			break;
		case GLUT_KEY_DOWN:
			beta -= 0.1f;
			if(beta < -1.5f)
				beta = -1.5;
			break;
		case GLUT_KEY_RIGHT:
			alpha -= 0.1f;
			break;
	}
	glutPostRedisplay();
}


//We assume that the .xml and .3d files passed are correct.
int main(int argc, char** argv){
	if(argc != 2){
		cerr << "Usage: " << argv[0] << " config_file\n";
		return 1;
	}

	TiXmlDocument doc(argv[1]);
	bool loadOkay = doc.LoadFile();
	if(!loadOkay){
		cout << "Error loading file '" << argv[1] << "'.\n";
		return 1;
	}

	string directory = directory_of_file(argv[1]);

	TiXmlHandle docHandle(&doc);
	TiXmlElement* model = doc.FirstChild("scene")->FirstChild("model")->ToElement();
	for(; model; model=model->NextSiblingElement()){
		string filename; 
		int r = model->QueryStringAttribute("file", &filename);
		
		if(r == TIXML_SUCCESS) {
			int n_vertex, n_triangles;
			ifstream file(directory + filename);
			figure triangles;

			if(!file) {
				cerr << "The file \"" << filename << "\" was not found.\n";
			}
			file >> n_vertex; // reads the number of vertices from the file
			n_triangles = n_vertex/3;
			
			for(int i = 0; i < n_triangles; i++){
				float color_r, color_g, color_b;
				point ps[3];
				for(int j = 0; j < 3; j++){
					float px, py, pz;
					file >> px;
					file >> py;
					file >> pz;
					point p(px, py, pz);
					ps[j] = p;
				}
				color_r = randFloat();
				color_g = randFloat();
				color_b = randFloat();

				triangle t(ps[0], ps[1], ps[2], color_r, color_g, color_b);
				triangles.push_back(t);
			}	
			file.close();
			figures.push_back(triangles);
		}
	}

	// init GLUT and the window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH|GLUT_DOUBLE|GLUT_RGBA);
	glutInitWindowPosition(100,100);
	glutInitWindowSize(800,800);
	glutCreateWindow("Pratical Assignment 1");
		
	// Required callback registry 
	glutDisplayFunc(renderScene);
	glutReshapeFunc(changeSize);
	
	// Callback registration for keyboard processing
	glutKeyboardFunc(processKeys);
	glutSpecialFunc(processSpecialKeys);

	//  OpenGL settings
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	
	// enter GLUT's main cycle
	glutMainLoop();

	return 1;
}
