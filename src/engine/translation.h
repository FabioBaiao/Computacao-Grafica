#ifndef TRANSLATION
#define TRANSLATION
#include "geoTransform.h"
#define TRANSLATE_UNIT 0.2

class translationCoords : public geoTransform {
public:
    float x, y, z;
    translationCoords(float fx, float fy, float fz) : x{fx}, y{fy}, z{fz} {};
    void apply() {
        glTranslatef(x, y, z);
    }
    void increaseX() {x += TRANSLATE_UNIT;}
    void increaseY() {y += TRANSLATE_UNIT;}
    void increaseZ() {z += TRANSLATE_UNIT;}
    void decreaseX() {x -= TRANSLATE_UNIT;}
    void decreaseY() {y -= TRANSLATE_UNIT;}
    void decreaseZ() {z -= TRANSLATE_UNIT;}
};

class translationTime : public geoTransform {
public:
    float time;
    float **points;
    int n;
    float up[3] = {0, 1, 0};
    translationTime(float ftime, float **fpoints, int in) : time{ftime*(float)1e3}, points{fpoints}, n{in} {};


    void renderCatmullRomCurve() {
        float pos[3];
        float deriv[3];

        // draw the curves using line segments - GL_LINE_LOOP
        // note: ONLY the models are required to be drawn with VBO's
        glBegin(GL_LINE_LOOP);
        //glColor3f(1,1,1);
        float spec[4] = {1.0f, 1.0f, 1.0f, 1.0f};
        glMaterialfv(GL_FRONT, GL_SPECULAR, spec);
        for(float i = 0.0; i < 1.0; i+=0.0001) {
            getGlobalCatmullRomPoint(i, (float *)pos, (float *)deriv);
            glVertex3f(pos[0], pos[1], pos[2]);
        }
        glEnd();
    }

    void multMatrixVector(float *m, float *v, float *res) {

        for (int j = 0; j < 4; ++j) {
            res[j] = 0;
            for (int k = 0; k < 4; ++k) {
                res[j] += v[k] * m[j * 4 + k];
            }
        }

    }

    void cross(float *a, float *b, float *res) {

        res[0] = a[1]*b[2] - a[2]*b[1];
        res[1] = a[2]*b[0] - a[0]*b[2];
        res[2] = a[0]*b[1] - a[1]*b[0];
    }


    void normalize(float *a) {

        float l = sqrt(a[0]*a[0] + a[1] * a[1] + a[2] * a[2]);
        a[0] = a[0]/l;
        a[1] = a[1]/l;
        a[2] = a[2]/l;
    }

    void getCatmullRomPoint(float t, float *p0, float *p1, float *p2, float *p3, float *pos, float *deriv) {

        // catmull-rom matrix
        float m[4][4] = {	{-0.5f,  1.5f, -1.5f,  0.5f},
            { 1.0f, -2.5f,  2.0f, -0.5f},
            {-0.5f,  0.0f,  0.5f,  0.0f},
            { 0.0f,  1.0f,  0.0f,  0.0f}
        };

        float px[4] = {p0[0], p1[0], p2[0], p3[0]};
        float py[4] = {p0[1], p1[1], p2[1], p3[1]};
        float pz[4] = {p0[2], p1[2], p2[2], p3[2]};

        // Compute A = M * P
        float ax[4];
        float ay[4];
        float az[4];

        multMatrixVector((float *)m, (float *)px, (float *)ax);
        multMatrixVector((float *)m, (float *)py, (float *)ay);
        multMatrixVector((float *)m, (float *)pz, (float *)az);

        // Compute pos = T * A
        float T[4] = {t*t*t, t*t, t, 1};

        pos[0] = T[0]*ax[0] + T[1]*ax[1] + T[2]*ax[2] + T[3]*ax[3];
        pos[1] = T[0]*ay[0] + T[1]*ay[1] + T[2]*ay[2] + T[3]*ay[3];
        pos[2] = T[0]*az[0] + T[1]*az[1] + T[2]*az[2] + T[3]*az[3];

        // compute deriv = T' * A
        float TT[4] = {3*(t*t), 2*t, 1, 0};

        deriv[0] = TT[0]*ax[0] + TT[1]*ax[1] + TT[2]*ax[2] + TT[3]*ax[3];
        deriv[1] = TT[0]*ay[0] + TT[1]*ay[1] + TT[2]*ay[2] + TT[3]*ay[3];
        deriv[2] = TT[0]*az[0] + TT[1]*az[1] + TT[2]*az[2] + TT[3]*az[3];
        // ...
    }

    void getGlobalCatmullRomPoint(float gt, float *pos, float *deriv) {

        float t = gt * n; // this is the real global t
        int index = floor(t);  // which segment
        t = t - index; // where within  the segment

        // indices store the points
        int indices[4];
        indices[0] = (index + n-1)%n;
        indices[1] = (indices[0]+1)%n;
        indices[2] = (indices[1]+1)%n;
        indices[3] = (indices[2]+1)%n;

        getCatmullRomPoint(t, points[indices[0]], points[indices[1]], points[indices[2]], points[indices[3]], pos, deriv);
    }

    void renderCatmullRom () {
        float pos[3];
        float deriv[3];

        float t = glutGet(GLUT_ELAPSED_TIME)/time;

        getGlobalCatmullRomPoint(t, (float *)pos, (float *)deriv);
        normalize((float *)deriv);

        float z[3];
        cross((float *)deriv, (float *)up, (float *)z);
        normalize((float *)z);

        cross((float *)z, (float *)deriv, (float *) up);
        normalize((float *)up);

        float m[4][4] = {{deriv[0], deriv[1], deriv[2], 0},
            {up[0], up[1], up[2], 0},
            {z[0], z[1], z[2], 0},
            {pos[0], pos[1], pos[2], 1}
        };

        glMultMatrixf((float *)m);
    }

    void apply() {
        renderCatmullRomCurve();
        renderCatmullRom();
    }
};

#endif
