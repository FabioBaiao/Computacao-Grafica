#define TIXML_USE_STL 
#define BUFF_SIZE 128
#include "tinyxml/tinyxml.h"
#include "point.h"
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

// Camera control
float r = 10;
float alpha;
float beta;

using namespace std;

// Structure to save figures to draw
vector< vector<point> > figures;

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

void renderScene(void) {

	// clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// set the camera
	glLoadIdentity();
	gluLookAt(r*cos(beta)*cos(alpha), r*sin(beta), r*cos(beta)*sin(alpha), 
		      0.0,0.0,0.0,
			  0.0f,1.0f,0.0f);

	//loop aqui
	for(auto figure:figures){
		glBegin(GL_TRIANGLES);
		// change color randomly?
		glColor3f(0,0,0);
		for(auto point:figure){
			glVertex3f(point.x, point.y, point.z);
		}
		glEnd();
	}
	// End of frame
	glutSwapBuffers();
}


void processKeys(unsigned char c, int xx, int yy) {
// put code to process regular keys in here
	switch(toupper(c)){
		case 'M': r+=0.2;
			  break;
		case 'L': r-=0.2;
			  break;
	}
	glutPostRedisplay();
}


void processSpecialKeys(int key, int xx, int yy) {
// put code to process special keys in here
	switch(key){
		case GLUT_KEY_UP:
			if(beta < M_PI/4-0.2)
				beta += 0.2;
			break;
		case GLUT_KEY_LEFT:
			alpha += 0.2;
			break;
		case GLUT_KEY_DOWN:
			if(beta > -M_PI/4+0.2)
				beta -= 0.2;
			break;
		case GLUT_KEY_RIGHT:
			alpha -= 0.2;
			break;
	}
	glutPostRedisplay();

}


//We assume that the .xml and .3d files passed are correct.
int main(int argc, char** argv){
	if(argc != 2){
		cout << "The command invocation must be in the form 'engine filename'\n";
		return -1;
	}

	TiXmlDocument doc(argv[1]);
	bool loadOkay = doc.LoadFile();
	if(!loadOkay){
		cout << "Error to load file " << argv[1] << ".\n";
		return -1;
	}

	// vector< vector<point> > figures1;
	// figures = figures1;

	TiXmlHandle docHandle(&doc);
	TiXmlElement* model = doc.FirstChild("scene")->FirstChild("model")->ToElement();
	for(; model; model=model->NextSiblingElement()){
		string filename; 
		int r = model->QueryStringAttribute("file", &filename);
		if(r == TIXML_SUCCESS){
			int n_vertex;
			ifstream file;
			file.open(filename);
			file >> n_vertex; // reads the number of vertices in a file
			std::vector<point> points;
			for(int i = 0; i < n_vertex; i++){
				float px, py, pz; 
				file >> px;
				file >> py;
				file >> pz;
				point p(px, py, px);
				points.push_back(p);
			}	
			file.close();
			figures.push_back(points);
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
