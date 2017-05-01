#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <utility>

#include <ctype.h>
#include <math.h>
#include <GL/glew.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "geoTransform.h"
#include "group.h"
#include "rotation.h"
#include "scale.h"
#include "tinyxml2.h"
#include "translation.h"

#define ANG2RAD M_PI/180

using namespace std;
using namespace tinyxml2;

// Scene models
typedef vector<float> Model;
vector<group> groups;

map<string, Model> models;
// maps model_name to buffer id and number of vertices
map<string, pair<GLuint,int>> model_to_buffer;

// maps keybindings to increase and decrease angle actions in rotations
map<char, rotation*> increaseBindings, decreaseBindings;
//vector<char> keysInUse{

// VBO's
int n_models;
GLuint* buffers;

// Camera control
float r = 10.0f;
float alpha;
float beta;

float pitch = 0.0f, yaw = 0.0f;

float Px = 0.0f, Py = 0.0f, Pz = 0.0f;
float lookX = Px + cos(pitch) * sin(yaw);
float lookY = Py + sin(pitch);
float lookZ = Pz + cos(pitch) * cos(yaw);

// Polygon Mode
GLenum modes[] = {GL_FILL, GL_LINE, GL_POINT};
GLenum mode;

// directory of the read file
string directory;

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
    gluPerspective(45.0f,ratio, 1.0f,1000.0f);
    glMatrixMode(GL_MODELVIEW);
}

void drawModel(string model, color c) {
    auto buffer_id_size = model_to_buffer[model];

    glColor3f(c.r, c.g, c.b);
    glBindBuffer(GL_ARRAY_BUFFER, buffer_id_size.first);
    glVertexPointer(3,GL_FLOAT, 0, 0);
    glDrawArrays(GL_TRIANGLES, 0, buffer_id_size.second/3);
}

void drawRandom(randomModel rnd) {
    for (auto spec : rnd.specs) {
        int n = spec.n;
        float maxR = spec.maxR;
        float minR = spec.minR;
        float maxS = spec.maxS;
        float minS = spec.minS;
        int i = 0;
        for (string model : rnd.models) {
            color c = rnd.modelsColor[i];
            for (int j = 0; j < n; j++) {
                float r = rand();
                float alfa = rand();
                float s = rand();
                // polar coordinates
                r = r/RAND_MAX * (maxR - minR) + minR;
                alfa = alfa/RAND_MAX * 2*M_PI;
                // random scale
                s = s/RAND_MAX * (maxS - minS) + minS;

                glPushMatrix();
                glTranslatef(r * sin(alfa), 0, r * cos(alfa));
                glScalef(s, s, s);
                drawModel(model, c);
                glPopMatrix();
            }
            i++;
        }
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

        drawModel(model, c);
        i++;
    }
    for (auto rnd : g.randoms) {
        drawRandom(rnd);
    }
    for(auto gr : g.childGroups) {
        drawGroup(gr);
    }
    glPopMatrix();
}

void renderScene(void) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPolygonMode(GL_FRONT, modes[mode]);
    glLoadIdentity();
    gluLookAt(Px, Py, Pz,
              lookX,lookY,lookZ,
              0.0f,1.0f,0.0f);

    srand(0);
    for(auto g : groups) {
        drawGroup(g);
    }
    glutSwapBuffers();
}

void lookAt(float alpha, float beta) {
    lookX = Px + cos(beta) * sin(alpha);
    lookY = Py + sin(beta);
    lookZ = Pz + cos(beta) * cos(alpha);
}

float magnitude(float vx, float vy, float vz) {
    return sqrtf((vx*vx) + (vy*vy) + (vz*vz));
}

void normalize(float* vx, float* vy, float* vz) {
    float mag = magnitude(*vx, *vy, *vz);
    *vx = (*vx) / mag;
    *vy = (*vy) / mag;
    *vz = (*vz) / mag;
}

void processKeys(unsigned char c, int xx, int yy) {
    float k = 0.5f;
    float upX = 0.0f, upY = 1.0f, upZ = 0.0f;

    /* forward vector */
    float dX = lookX - Px;
    float dY = lookY - Py;
    float dZ = lookZ - Pz;

    normalize(&dX, &dY, &dZ);
    float rX, rY, rZ;
    char cc = toupper(c);
    switch(cc) {
    case 27:
        /* ESC key */
        exit(0);
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
        /* cross product: forwardV x up */
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
    default:
        auto itDec = decreaseBindings.find(cc);

        if(itDec != decreaseBindings.end()) {
            itDec->second->decreaseAngle();
        } else{ 
            auto itInc = increaseBindings.find(cc);
            if(itInc != increaseBindings.end()) {
                 itInc->second->increaseAngle();
            }
        }
        break;
    }
    glutPostRedisplay();
}

void processSpecialKeys(int key, int xx, int yy) {
    switch(key) {
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

void parsePoint(float **points, int i, XMLElement* elem) {
    float x, y, z;
    elem->QueryFloatAttribute("X", &x);
    elem->QueryFloatAttribute("Y", &y);
    elem->QueryFloatAttribute("Z", &z);

    points[i] = (float *) malloc(sizeof(float) * 3);

    points[i][0] = x;
    points[i][1] = y;
    points[i][2] = z;
}

void parseTranslate(group& g, XMLElement* elem) {
    float x, y, z, time;
    x = y = z = time = 0.0f;

    XMLError hasTime = elem->QueryFloatAttribute("time", &time);
    
    if(hasTime != XML_NO_ATTRIBUTE) {
        // animated translation
        XMLElement *child = elem->FirstChildElement();
        float **points = NULL;
        int n = 0;
        for(; child; child = child->NextSiblingElement()) {
            n++;
            points = (float **) realloc(points, n*sizeof(float *));
            parsePoint(points, n-1, child);
        }
        g.transforms.push_back(new translationTime(time, points, n));
    } else {
        // static translation
    	elem->QueryFloatAttribute("X", &x);
    	elem->QueryFloatAttribute("Y", &y);
    	elem->QueryFloatAttribute("Z", &z);
    	g.transforms.push_back(new translationCoords(x,y,z));
    }
}

void parseRotate(group& g, XMLElement* elem) {
    float angle, time, axisX, axisY, axisZ;
    angle = axisX = axisY = axisZ = time = 0.0f;

    elem->QueryFloatAttribute("time", &time);
    elem->QueryFloatAttribute("angle", &angle);
    elem->QueryFloatAttribute("axisX", &axisX);
    elem->QueryFloatAttribute("axisY", &axisY);
    elem->QueryFloatAttribute("axisZ", &axisZ);


    rotation* rot = new rotation(angle, axisX, axisY, axisZ, time);
    g.transforms.push_back(rot);

    const char * decreaseBind = elem->Attribute("decreaseBind");
    const char * increaseBind = elem->Attribute("increaseBind");

    if(decreaseBind && strlen(decreaseBind) == 1) {
        // if there is a keybinding to lower the angle
        char c = toupper(decreaseBind[0]);
        decreaseBindings[c] = rot;
    } 

    if(increaseBind && strlen(increaseBind) == 1) {
        // if there is a keybinding to lower the angle
        char c = toupper(increaseBind[0]);
        increaseBindings[c] = rot;
    } 

}

void parseScale(group& g, XMLElement* elem) {
    float x, y, z;
    x = y = z = 1.0f;
    elem->QueryFloatAttribute("X", &x);
    elem->QueryFloatAttribute("Y", &y);
    elem->QueryFloatAttribute("Z", &z);
    g.transforms.push_back(new scale(x,y,z));
}

color& parseColor(XMLElement* model) {
    float rr, gg, bb;
    int r_r, r_g, r_b;
    rr = gg = bb = 0.0;
    r_r = model->QueryFloatAttribute("R", &rr);
    r_g = model->QueryFloatAttribute("G", &gg);
    r_b = model->QueryFloatAttribute("B", &bb);
    if(r_r != XML_SUCCESS && r_g != XML_SUCCESS && r_b != XML_SUCCESS) {
        // the color defaults to white if not specified
        rr = gg = bb = 1;
    }
    color* c = new color(rr,gg,bb);
    return *c;
}

void readFile(ifstream& file, string fName) {
    int n_vertex, n_values;
    Model model_read;

    file >> n_vertex;
    n_values = 3*n_vertex;

    for(int i = 0; i < n_values; i++) {
        float tmp;
        file >> tmp;
        model_read.push_back(tmp);
    }
    file.close();
    models[fName] = model_read;
}

void parseModel(group& g, XMLElement * model) {
    const char* filename= model->Attribute("file");
    if(filename != NULL) {
        string fName = string(filename);
        ifstream file(directory + filename);
        if(!file) {
            cerr << "The file \"" << filename << "\" was not found.\n";
        }
        color c = parseColor(model);
        g.models.push_back(fName);
        g.modelsColor.push_back(c);
        if(models.find(fName) == models.end()) {
            // if the model file was not already read
            readFile(file, fName);
        }
    }
}

void parseModel(randomModel& r, XMLElement * model) {
    const char* filename= model->Attribute("file");
    if(filename != NULL) {
        string fName = string(filename);
        ifstream file(directory + filename);
        if(!file) {
            cerr << "The file \"" << filename << "\" was not found.\n";
        }
        color c = parseColor(model);
        r.models.push_back(fName);
        r.modelsColor.push_back(c);
        if(models.find(fName) == models.end()) {
            // if the model file was not already read
            readFile(file, fName);
        }
    }
}

void parseSpecs(randomModel& r, XMLElement * model) {
    int n = 0;
    float minR=0.0, maxR=0.0, minS=1.0, maxS=1.0;
    minR = maxR = minS = maxS = 0;
    model->QueryIntAttribute("N", &n);
    model->QueryFloatAttribute("minRadius", &minR);
    model->QueryFloatAttribute("maxRadius", &maxR);
    model->QueryFloatAttribute("minScale", &minS);
    model->QueryFloatAttribute("maxScale", &maxS);

    randSpecs* rs = new randSpecs(n, minR, maxR, minS, maxS);
    r.specs.push_back(*rs);
}

void parseRandom(group& g, XMLElement * elem) {
    randomModel r;
    XMLElement *child = elem->FirstChildElement();
    for(; child; child = child->NextSiblingElement()) {
        string type = string(child->Name());
        if (type == "model") {
            parseModel(r, child);
        }
        else if(type == "specs") {
            parseSpecs(r, child);
        }
    }
    g.randoms.push_back(r);
}

void parseModels(group& g, XMLElement * elem) {
    XMLElement *child = elem->FirstChildElement();
    for(; child; child = child->NextSiblingElement()) {
        string type = string(child->Name());
        if (type == "model") {
            parseModel(g, child);
        }
        else if(type == "random") {
            parseRandom(g, child);
        }
    }
}

group parseGroup(XMLElement *gr) {
    group g;
    XMLElement *child = gr->FirstChildElement();
    for( ; child; child = child->NextSiblingElement()) {
        string type = string(child->Name());
        if(type == "translate") {
            parseTranslate(g, child);
        } else if(type == "rotate") {
            parseRotate(g, child);
        } else if(type == "scale") {
            parseScale(g, child);
        } else if(type == "group") {
            group g_child = parseGroup(child);
            g.childGroups.push_back(g_child);
        } else if(type == "models") {
            parseModels(g, child);
        }
    }
    return g;
}

void parseInitialPosition(XMLElement* scene) {
    // if a coordinate is not specified, it will default to 0
    scene->QueryFloatAttribute("camX", &Px);
    scene->QueryFloatAttribute("camY", &Py);
    scene->QueryFloatAttribute("camZ", &Pz);
    scene->QueryFloatAttribute("yaw", &yaw);
    scene->QueryFloatAttribute("pitch", &pitch);
    yaw *= ANG2RAD;
    pitch *= ANG2RAD;
    lookX = Px + cos(pitch) * sin(yaw);
    lookY = Py + sin(pitch);
    lookZ = Pz + cos(pitch) * cos(yaw);
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

    parseInitialPosition(scene);

    XMLElement* gr = scene->FirstChildElement("group");

    /* Read all the groups */
    for(; gr; gr = gr->NextSiblingElement()) {
        group g = parseGroup(gr);
        groups.push_back(g);
    }

    // init GLUT and the window
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH|GLUT_DOUBLE|GLUT_RGBA);
    glutInitWindowPosition(100,100);
    glutInitWindowSize(glutGet(GLUT_SCREEN_WIDTH), glutGet(GLUT_SCREEN_HEIGHT));
    glutCreateWindow("Pratical Assignment");

    // Callback registry
    glutDisplayFunc(renderScene);
    glutIdleFunc(renderScene);
    glutReshapeFunc(changeSize);
    glutKeyboardFunc(processKeys);
    glutSpecialFunc(processSpecialKeys);

    glewInit();

    n_models = models.size();
    buffers = (GLuint *) malloc(sizeof(GLuint) * n_models);
    glGenBuffers(n_models, buffers);

    int i = 0;
    for(auto it = models.begin(); it != models.end(); it++) {
        auto model = it->second;

        /* maps the name of the model to the id of the buffer and number of vertices */
        pair<GLuint, int> id_number_elements(buffers[i], model.size());
        model_to_buffer[it->first] = id_number_elements;

        /* fills the buffers */
        glBindBuffer(GL_ARRAY_BUFFER, buffers[i++]);
        glBufferData(GL_ARRAY_BUFFER, model.size() * sizeof(float), model.data(), GL_STATIC_DRAW);
    }

    //  OpenGL settings
    glEnable(GL_DEPTH_TEST);
    glMatrixMode(GL_MODELVIEW);
    glEnable(GL_CULL_FACE);
    glEnableClientState(GL_VERTEX_ARRAY);

    glutMainLoop();

    return 1;
}
