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

/* Size (in bytes) of the buffer used to store each line of a .patch file (one at a time) */
#define BUFF_SIZE 1024

/* s coordinate of the center of the upper base in a standard texture of a frustum */
#define S_FRUSTUM_UBASE_CENTER 0.4375f
/* same as the previous macro, but for the lower base */
#define S_FRUSTUM_LBASE_CENTER 0.8125f
/* t coordinate of the center of each base in a standard texture of a frustum (both bases' centers have the same t)s */
#define T_FRUSTUM_BASES_CENTER 0.1875f
/* base radius of each base in a standard texture for a frustum */
#define FRUSTUM_TEX_BASE_RADIUS 0.1875f
/* lower t coordinate of the start of the side area in a standard texture of a frustum */
#define LOW_T_FRUSTUM_SIDE_AREA 0.375f

typedef struct point {
    float x;
    float y;
    float z;
} point;

point *cpoints;     // control points
int **indexes;      // indexes of the points for each patch
int patches;        // number of patches
int ncpoints;       // number of control points

void cross(float *a, float *b, float *res) {

    res[0] = a[1] * b[2] - a[2] * b[1];
    res[1] = a[2] * b[0] - a[0] * b[2];
    res[2] = a[0] * b[1] - a[1] * b[0];
}

float length(float *v) {
	float res = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    
    return res;
}

void normalize(float *a) {
    float l = length(a);
    if (l == 0.0f)
        return;
    
    a[0] = a[0] / l;
    a[1] = a[1] / l;
    a[2] = a[2] / l;
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
    for (int i = 0; i < 4; ++i) {
        res[i] = 0;
        for (int j = 0; j < 4; ++j) {
            res[i] += v[j] * m[j * 4 + i];
        }
    }
}

void multMatrixMatrix(float *m1, float *m2, float *res) {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            res[i * 4 + j] = 0.0f;
            for (int k = 0; k < 4; ++k)
                res[i * 4 + j] += m1[i * 4 + k] * m2[k * 4 + j];
        }
    }
}

string planeXZ (float x, float y, float z, int nDiv, int orient) {
    ostringstream os;
    int count = 0;
    float deltaX = x / nDiv;
    float deltaZ = z / nDiv;

    x = x / 2;
    z = z / 2;

    float xi = x;
    float zi = z;

    for (int i = 1; i <= nDiv; i++) {
        for (int j = 1; j <= nDiv; j++) {
            if (orient == UPWARDS) {
                os << xi << " " << y << " " << zi << '\n';
                os << xi << " " << y << " " << (zi - deltaZ) << '\n';
                os << (xi - deltaX) << " " << y << " " << (zi - deltaZ) << '\n';

                os << xi << " " << y << " " << zi << '\n';
                os << (xi - deltaX) << " " << y << " " << (zi - deltaZ) << '\n';
                os << (xi - deltaX) << " " << y << " " << zi << '\n';
            }
            else {
                os << xi << " " << y << " " << zi << '\n';
                os << (xi - deltaX) << " " << y << " " << zi << '\n';
                os << (xi - deltaX) << " " << y << " " << (zi - deltaZ) << '\n';

                os << xi << " " << y << " " << zi << '\n';
                os << (xi - deltaX) << " " << y << " " << (zi - deltaZ) << '\n';
                os << xi << " " << y << " " << (zi - deltaZ) << '\n';
            }
            count += 6;
            xi = x - j * deltaX;
        }
        zi = z - i * deltaZ;
        xi = x;
    }
    return to_string(count) + "\n" + os.str();
}

string planeXZnormals (int nDiv, int orient) {
    ostringstream os;
    for (int i = 1; i <= nDiv; i++) {
        for (int j = 1; j <= nDiv; j++) {
            for (int k = 1; k <= 6; k++) {
                os << "0.0 " << orient * 1.0f << " 0.0\n";
            }
        }
    }
    return os.str();
}

string planeXY (float x, float y, float z, int nDiv, int orient) {
    ostringstream os;
    int count = 0;
    float deltaX = x / nDiv;
    float deltaY = y / nDiv;

    x = x / 2;
    y = y / 2;

    float xi = x;
    float yi = -y;

    for (int i = 1; i <= nDiv; i++) {
        for (int j = 1; j <= nDiv; j++) {
            if (orient == UPWARDS) {
                os << xi << " " << yi << " " << z << '\n';
                os << xi << " " << (yi + deltaY) << " " << z << '\n';
                os << (xi - deltaX) << " " << (yi + deltaY) << " " << z << '\n';

                os << xi << " " << yi << " " << z << '\n';
                os << (xi - deltaX) << " " << (yi + deltaY) << " " << z << '\n';
                os << (xi - deltaX) << " " << yi << " " << z << '\n';
            }
            else {
                os << xi << " " << yi << " " << z << '\n';
                os << (xi - deltaX) << " " << yi << " " << z << '\n';
                os << (xi - deltaX) << " " << (yi + deltaY) << " " << z << '\n';

                os << xi << " " << yi << " " << z << '\n';
                os << (xi - deltaX) << " " << (yi + deltaY) << " " << z << '\n';
                os << xi << " " << (yi + deltaY) << " " << z << '\n';
            }
            count += 6;
            xi = x - j * deltaX;
        }
        yi = -y + i * deltaY;
        xi = x;
    }
    return to_string(count) + "\n" + os.str();
}

string planeXYnormals (int nDiv, int orient) {
    ostringstream os;
    for (int i = 1; i <= nDiv; i++) {
        for (int j = 1; j <= nDiv; j++) {
            for (int k = 1; k <= 6; k++) {
                os << "0.0 0.0 " << orient * 1.0f << "\n";
            }
        }
    }
    return os.str();
}

string planeYZ (float x, float y, float z, int nDiv, int orient) {
    ostringstream os;
    int count = 0;
    float deltaY = y / nDiv;
    float deltaZ = z / nDiv;

    y = y / 2;
    z = z / 2;

    float yi = -y;
    float zi = -z;

    for (int i = 1; i <= nDiv; i++) {
        for (int j = 1; j <= nDiv; j++) {
            if (orient == UPWARDS) {
                os << x << " " << yi << " " << zi << '\n';
                os << x << " " << (yi + deltaY) << " " << zi << '\n';
                os << x << " " << (yi + deltaY) << " " << (zi + deltaZ) << '\n';
                os << x << " " << yi << " " << zi << '\n';
                os << x << " " << (yi + deltaY) << " " << (zi + deltaZ) << '\n';
                os << x << " " << yi << " " << (zi + deltaZ) << '\n';
            } else {
                os << x << " " << yi << " " << zi << '\n';
                os << x << " " << yi << " " << (zi + deltaZ) << '\n';
                os << x << " " << (yi + deltaY) << " " << (zi + deltaZ) << '\n';
                os << x << " " << yi << " " << zi << '\n';
                os << x << " " << (yi + deltaY) << " " << (zi + deltaZ) << '\n';
                os << x << " " << (yi + deltaY) << " " << zi << '\n';
            }
            count += 6;
            zi = -z + j * deltaZ;
        }
        yi = -y + i * deltaY;
        zi = -z;
    }
    return to_string(count) + "\n" + os.str();
}

string planeYZnormals (int nDiv, int orient) {
    ostringstream os;
    for (int i = 1; i <= nDiv; i++) {
        for (int j = 1; j <= nDiv; j++) {
            for (int k = 1; k <= 6; k++) {
                os << orient * 1.0f << " 0.0 0.0\n";
            }
        }
    }
    return os.str();
}

string planeTexCoords(int nDiv, int orient) {
    ostringstream os;
    float deltaS = 1.0f / nDiv;
    float deltaT = 1.0f / nDiv;

    float s = 1.0f;
    float t = 0.0f;

    for (int i = 1; i <= nDiv; i++) {
        for (int j = 1; j <= nDiv; j++) {
            if (orient == UPWARDS) {
                os << s << " " << t << "\n";
                os << s << " " << (t + deltaT) << "\n";
                os << (s - deltaS) << " " << (t + deltaT) << "\n";

                os << s << " " << t << "\n";
                os << (s - deltaS) << " " << (t + deltaT) << "\n";
                os << (s - deltaS) << " " << t << "\n";
            }
            else {
                os << s << " " << t << '\n';
                os << (s - deltaS) << " " << t << '\n';
                os << (s - deltaS) << " " << (t + deltaT) << '\n';

                os << s << " " << t << '\n';
                os << (s - deltaS) << " " << (t + deltaT) << '\n';
                os << s << " " << (t + deltaT) << '\n';
            }
            s = 1.0f - j * deltaS;
        }
        t = i * deltaT;
        s = 1.0f;
    }
    return os.str();
}

string boxTexCoords(int nDiv) {
    string upwardsTexCoords = planeTexCoords(nDiv, UPWARDS);
    string downwardsTexCoords = planeTexCoords(nDiv, DOWNWARDS);

    return upwardsTexCoords + downwardsTexCoords + upwardsTexCoords +
           downwardsTexCoords + upwardsTexCoords + downwardsTexCoords;
}

string box(float x, float y, float z, int nDiv) {
    int count;
    string plane, line;
    ostringstream os;
    istringstream iss;

    plane = planeXY(x, y, z / 2, nDiv, UPWARDS);
    istringstream iss1 (plane);
    getline(iss1, line);
    count = stoi(line, nullptr, 10);
    while (getline(iss1, line)) {
        os << line << '\n';
    }
    plane = planeXY(x, y, -z / 2, nDiv, DOWNWARDS);
    istringstream iss2 (plane);
    getline(iss2, line);
    count += stoi(line, nullptr, 10);
    while (getline(iss2, line)) {
        os << line << '\n';
    }
    plane = planeXZ(x, y / 2, z, nDiv, UPWARDS);
    istringstream iss3 (plane);
    getline(iss3, line);
    count += stoi(line, nullptr, 10);
    while (getline(iss3, line)) {
        os << line << '\n';
    }
    plane = planeXZ(x, -y / 2, z, nDiv, DOWNWARDS);
    istringstream iss4 (plane);
    getline(iss4, line);
    count += stoi(line, nullptr, 10);
    while (getline(iss4, line)) {
        os << line << '\n';
    }
    plane = planeYZ(x / 2, y, z, nDiv, UPWARDS);
    istringstream iss5 (plane);
    getline(iss5, line);
    count += stoi(line, nullptr, 10);
    while (getline(iss5, line)) {
        os << line << '\n';
    }
    plane = planeYZ(-x / 2, y, z, nDiv, DOWNWARDS);
    istringstream iss6 (plane);
    getline(iss6, line);
    count += stoi(line, nullptr, 10);
    while (getline(iss6, line)) {
        os << line << '\n';
    }
    return to_string(count) + "\n" + os.str();
}

string boxNormals(int nDiv) {
    string normals, line;
    ostringstream os;
    istringstream iss;

    normals = planeXYnormals(nDiv, UPWARDS);
    istringstream iss1(normals);
    while (getline(iss1, line)) {
        os << line << '\n';
    }
    normals = planeXYnormals(nDiv, DOWNWARDS);
    istringstream iss2 (normals);
    while (getline(iss2, line)) {
        os << line << '\n';
    }
    normals = planeXZnormals(nDiv, UPWARDS);
    istringstream iss3 (normals);
    while (getline(iss3, line)) {
        os << line << '\n';
    }
    normals = planeXZnormals(nDiv, DOWNWARDS);
    istringstream iss4 (normals);
    while (getline(iss4, line)) {
        os << line << '\n';
    }
    normals = planeYZnormals(nDiv, UPWARDS);
    istringstream iss5 (normals);
    while (getline(iss5, line)) {
        os << line << '\n';
    }
    normals = planeYZnormals(nDiv, DOWNWARDS);
    istringstream iss6 (normals);
    while (getline(iss6, line)) {
        os << line << '\n';
    }
    return os.str();
}

string annulus(float dist, float smj, float smn, int slices) {
    int nPoints = 0;
    float deltaAlpha;
    ostringstream os;

    deltaAlpha = (2.0f * M_PI) / slices;
    for (int i = 0; i < slices; ++i) {
        float alpha = i * deltaAlpha;
        float nextAlpha = alpha + deltaAlpha;

        os << smj * cosf(alpha) << " 0.0 " << smn * sinf(alpha) << '\n';
        os << (smj + dist) * cosf(nextAlpha) << " 0.0 " << (smn + dist) * sinf(nextAlpha) << '\n';
        os << (smj + dist) * cosf(alpha) << " 0.0 " << (smn + dist) * sinf(alpha) << '\n';
        os << smj * cosf(alpha) << " 0.0 " << smn * sinf(alpha) << '\n';
        os << smj * cosf(nextAlpha) << " 0.0 " << smn * sinf(nextAlpha) << '\n';
        os << (smj + dist) * cosf(nextAlpha) << " 0.0 " << (smn + dist) * sinf(nextAlpha) << '\n';
        nPoints += 6;
    }
    return (to_string(nPoints) + "\n" + os.str());
}

string annulusNormals(int slices) {
    ostringstream os;
    for (int i = 0; i < slices; ++i) {
        for (int j = 0; j < 6; ++j) {
            os << "0.0 1.0 0.0\n";
        }
    }
    return os.str();
}

string annulusTexCoords(int slices) {
    ostringstream os;

    for (int i = 0; i < slices; ++i) {
        os << "1.0 0.0\n";
        os << "0.0 1.0\n";
        os << "0.0 0.0\n";
        os << "1.0 0.0\n";
        os << "1.0 1.0\n";
        os << "0.0 1.0\n";
    }
    return os.str();
}

string ellipsoid(float a, float b, float c, int stacks, int slices) {
    int nPoints = 0;
    float deltaAlpha = 2.0f * M_PI / slices;
    float deltaBeta = M_PI / stacks;
    ostringstream os;

    for (int i = 0; i < stacks; i++) {
        float beta = i * deltaBeta;
        float nextBeta = beta + deltaBeta;

        for (int j = 0; j < slices; j++) {
            float alpha = j * deltaAlpha;
            float nextAlpha = alpha + deltaAlpha;

            if (i < stacks - 1) {
                os << a * sinf(beta) * sinf(alpha) << " " << b * cosf(beta) << " " << c * sinf(beta) * cosf(alpha) << '\n';
                os << a * sinf(nextBeta) * sinf(alpha) << " " << b * cosf(nextBeta) << " " << c * sinf(nextBeta) * cosf(alpha) << '\n';
                os << a * sinf(nextBeta) * sinf(nextAlpha) << " " << b * cosf(nextBeta) << " " << c * sinf(nextBeta) * cosf(nextAlpha) << '\n';
                nPoints += 3;
            }
            if (i > 0) {
                os << a * sinf(beta) * sinf(alpha) << " " << b * cosf(beta) << " " << c * sinf(beta) * cosf(alpha) << '\n';
                os << a * sinf(nextBeta) * sinf(nextAlpha) << " " << b * cosf(nextBeta) << " " << c * sinf(nextBeta) * cosf(nextAlpha) << '\n';
                os << a * sinf(beta) * sinf(nextAlpha) << " " << b * cosf(beta) << " " << c * sinf(beta) * cosf(nextAlpha) << '\n';
                nPoints += 3;
            }
        }

    }
    return (to_string(nPoints) + "\n" + os.str());
}

/** 
 * Function that given a, b and c of an ellipsoid, alpha and beta of a surface point and an ostringstream,
 * appends the coordinates of the normal of that surface point to the provided ostringstream.
 */
void appendEllipsoidNormal(float a, float b, float c, float alpha, float beta, ostringstream& os) {
	float n[3] = {
		sinf(beta) * sinf(alpha) / a,
		cosf(beta) / b,
		sinf(beta) * cosf(alpha) / c
	};
	normalize(n);

	os << n[0] << ' ' << n[1] << ' ' << n[2] << '\n';
}

string ellipsoidNormals(float a, float b, float c, int stacks, int slices) {
    float deltaAlpha = 2.0f * M_PI / slices;
    float deltaBeta = M_PI / stacks;
    ostringstream os;

    for (int i = 0; i < stacks; i++) {
        float beta = i * deltaBeta;
        float nextBeta = beta + deltaBeta;

        for (int j = 0; j < slices; j++) {
            float alpha = j * deltaAlpha;
            float nextAlpha = alpha + deltaAlpha;

            if (i < stacks - 1) {
            	appendEllipsoidNormal(a, b, c, alpha, beta, os);
            	appendEllipsoidNormal(a, b, c, alpha, nextBeta, os);
            	appendEllipsoidNormal(a, b, c, nextAlpha, nextBeta, os);
            }
            if (i > 0) {
                appendEllipsoidNormal(a, b, c, alpha, beta, os);
                appendEllipsoidNormal(a, b, c, nextAlpha, nextBeta, os);
                appendEllipsoidNormal(a, b, c, nextAlpha, beta, os);
            }
        }
    }
    return os.str();
}

string ellipsoidTexCoords(int stacks, int slices) {
    ostringstream os;
    float deltaS = 1.0f / slices;
    float deltaT = 1.0f / stacks;

    for (int i = 0; i < stacks; ++i) {
        float t = (stacks - i) * deltaT;

        for (int j = 0; j < slices; ++j) {
            float s = j * deltaS;

            if (i < stacks - 1) {
                os << s << ' ' << t << '\n';
                os << s << ' ' << (t - deltaT)  << '\n';
                os << (s + deltaS) << ' ' << (t - deltaT) << '\n';
            }
            if (i > 0) {
                os << s << ' ' << t << '\n';
                os << (s + deltaS) << ' ' << (t - deltaT) << '\n';
                os << (s + deltaS) << ' ' << t << '\n';
            }
        }
    }
    return os.str();
}

string frustum(float baseRadius, float topRadius, float height, int slices, int stacks) {
    int i, j, nPoints = 0;
    float *cosCache, *sinCache;
    float dr, dh, alpha, dAlpha;
    ostringstream os;

    cosCache = new float[slices + 1];
    sinCache = new float[slices + 1];
    alpha = 0.0f;
    dAlpha = 2.0f * M_PI / slices;
    for (i = 0; i < slices; i++) {
        alpha = i * dAlpha;
        cosCache[i] = cosf(alpha);
        sinCache[i] = sinf(alpha);
    }
    cosCache[slices] = cosCache[0]; // cos(2 * M_PI) = cos(0)
    sinCache[slices] = sinCache[0]; // sin(2 * M_PI) = sin(0)

    if (baseRadius > 0.0f) {
        for (i = 0; i < slices; i++) { // lower base
            os << (baseRadius * sinCache[i]) << " 0.0 " << (baseRadius * cosCache[i]) << '\n';
            os << "0.0 0.0 0.0\n";
            os << (baseRadius * sinCache[i + 1]) << " 0.0 " << (baseRadius * cosCache[i + 1]) << '\n';
            nPoints += 3;
        }
    }
    if (topRadius > 0.0f) {
        for (i = 0; i < slices; ++i) { // upper base
            os << (topRadius * sinCache[i]) << ' ' << height << ' ' << (topRadius * cosCache[i]) << '\n';
            os << (topRadius * sinCache[i + 1]) << ' ' << height << ' ' << (topRadius * cosCache[i + 1]) << '\n';
            os << "0.0 " << height << " 0.0\n";
            nPoints += 3;
        }
    }
    dr = (baseRadius - topRadius) / stacks;
    dh = height / stacks;
    // side of the frustum
    for (i = 0; i < stacks; ++i) {
        float rLow, rHigh;
        float yLow, yHigh;

        rLow = baseRadius - i * dr;
        rHigh = rLow - dr;
        yLow = i * dh;
        yHigh = yLow + dh;
        if (rHigh > 0.0f) {
            for (j = 0; j < slices; ++j) {
                os << (rLow * sinCache[j]) << ' ' << yLow << ' ' << (rLow * cosCache[j]) << '\n';
                os << (rHigh * sinCache[j + 1]) << ' ' << yHigh << ' ' << (rHigh * cosCache[j + 1]) << '\n';
                os << (rHigh * sinCache[j]) << ' ' << yHigh << ' ' << (rHigh * cosCache[j]) << '\n';
                nPoints += 3;
            }
        }
        if (rLow > 0.0f) {
            for (j = 0; j < slices; ++j) {
                os << (rLow * sinCache[j]) << ' ' << yLow << ' ' << (rLow * cosCache[j]) << '\n';
                os << (rLow * sinCache[j + 1]) << ' ' << yLow << ' ' << (rLow * cosCache[j + 1]) << '\n';
                os << (rHigh * sinCache[j + 1]) << ' ' << yHigh << ' ' << (rHigh * cosCache[j + 1]) << '\n';
                nPoints += 3;
            }
        }
    }
    delete[] cosCache; delete[] sinCache;

    return (to_string(nPoints) + "\n" + os.str());
}

string frustumNormals(float baseRadius, float topRadius, float height, int slices, int stacks) {
    int i, j;
    float nxz, ny;
    float *cosCache, *sinCache;
    float dr, alpha, dAlpha;
    ostringstream os;

    cosCache = new float[slices + 1];
    sinCache = new float[slices + 1];
    alpha = 0.0f;
    dAlpha = 2.0f * M_PI / slices;
    for (i = 0; i < slices; ++i) {
        alpha = i * dAlpha;
        cosCache[i] = cosf(alpha);
        sinCache[i] = sinf(alpha);
    }
    cosCache[slices] = cosCache[0]; // cos(2 * M_PI) = cos(0)
    sinCache[slices] = sinCache[0]; // sin(2 * M_PI) = sin(0)

    if (baseRadius > 0.0f) {
        for (i = 0; i < slices; ++i) { // lower base
            for (j = 0; j < 3; ++j) {
                os << "0.0 -1.0 0.0\n";
            }
        }
    }
    if (topRadius > 0.0f) {
        for (i = 0; i < slices; ++i) { // upper base
            for (j = 0; j < 3; ++j) {
                os << "0.0 1.0 0.0\n";
            }
        }
    }
    
    if (baseRadius > topRadius) {
        float appexHeight = height / (1 - topRadius / baseRadius);
        float hypot = sqrt(baseRadius * baseRadius + appexHeight * appexHeight); // length of the hypotnuse
        
        ny = baseRadius / hypot; // y component of side area normal
        nxz = appexHeight / hypot; // sqrt(nx * nx + nz * nz) for each side area normal
    } else if (baseRadius < topRadius) {
        float appexHeight = height / (1 - baseRadius / topRadius);
        float hypot = sqrt(baseRadius * baseRadius + appexHeight * appexHeight); // length of the hypotnuse

        ny = -baseRadius / hypot; // y component of side area normal points down
        nxz = appexHeight / hypot;
    } else { // cylinder
        ny = 0.0f;
        nxz = 1.0f;
    }
    dr = (baseRadius - topRadius) / stacks;
    // side of the frustum
    for (i = 0; i < stacks; ++i) {
        float rLow, rHigh;

        rLow = baseRadius - i * dr;
        rHigh = rLow - dr;
        if (rHigh > 0.0f) {
            for (j = 0; j < slices; ++j) {
                os << (nxz * sinCache[j]) << ' ' << ny << ' ' << (nxz * cosCache[j]) << '\n';
                os << (nxz * sinCache[j + 1]) << ' ' << ny << ' ' << (nxz * cosCache[j + 1]) << '\n';
                os << (nxz * sinCache[j]) << ' ' << ny << ' ' << (nxz * cosCache[j]) << '\n';
            }
        }
        if (rLow > 0.0f) {
            for (j = 0; j < slices; ++j) {
                os << (nxz * sinCache[j]) << ' ' << ny << ' ' << (nxz * cosCache[j]) << '\n';
                os << (nxz * sinCache[j + 1]) << ' ' << ny << ' ' << (nxz * cosCache[j + 1]) << '\n';
                os << (nxz * sinCache[j + 1]) << ' ' << ny << ' ' << (nxz * cosCache[j + 1]) << '\n';
            }
        }
    }
    delete[] cosCache; delete[] sinCache;

    return os.str();
}

string frustumTexCoords(float baseRadius, float topRadius, int slices, int stacks) {
    int i, j;
    float alpha, deltaAlpha;
    float sCenter, tCenter, rTex;
    ostringstream os;

    deltaAlpha = 2.0f * M_PI / slices;
    rTex = FRUSTUM_TEX_BASE_RADIUS;
    tCenter = T_FRUSTUM_BASES_CENTER; // t coordinate is the same for both bases

    if (baseRadius > 0.0f) { // lower base texture coordinates
        sCenter = S_FRUSTUM_LBASE_CENTER;
        for (i = 0; i < slices; ++i) {
            alpha = i * deltaAlpha;
            os << (sCenter + rTex * sinf(alpha)) << ' ' << (tCenter - rTex * cosf(alpha)) << '\n';
            os << sCenter << ' ' << tCenter << '\n';
            os << (sCenter + rTex * sinf(alpha + deltaAlpha)) << ' ' << (tCenter - rTex * cosf(alpha + deltaAlpha)) << '\n';
        }
    }
    if (topRadius > 0.0f) {
        sCenter = S_FRUSTUM_UBASE_CENTER;
        for (i = 0; i < slices; ++i) { // upper base texture coordinates
            alpha = i * deltaAlpha;
            os << (sCenter + rTex * sinf(alpha)) << ' ' << (tCenter - rTex * cosf(alpha)) << '\n';
            os << (sCenter + rTex * sinf(alpha + deltaAlpha)) << ' ' << (tCenter - rTex * cosf(alpha + deltaAlpha)) << '\n';
            os << sCenter << ' ' << tCenter << '\n';
        }
    }

    float dr = (baseRadius - topRadius) / stacks;
    float deltaS = 1.0f / slices;
    float deltaT = (1.0f - LOW_T_FRUSTUM_SIDE_AREA) / stacks;
    for (i = 0; i < stacks; ++i) {
        float rLow = baseRadius - i * dr;
        float rHigh = rLow - dr;
        float s, t = LOW_T_FRUSTUM_SIDE_AREA + i * deltaT;

        if (rHigh > 0.0f) {
            for (j = 0; j < slices; ++j) {
                s = j * deltaS;
                os << s << ' ' << t << '\n';
                os << (s + deltaS) << ' ' << (t + deltaT) << '\n';
                os << s << ' ' << (t + deltaT) << '\n';
            }
        }
        if (rLow > 0.0f) {
            for (j = 0; j < slices; ++j) {
                s = j * deltaS;
                os << s << ' ' << t << '\n';
                os << (s + deltaS) << ' ' << t << '\n';
                os << (s + deltaS) << ' ' << (t + deltaT) << '\n';
            }
        }
    }
    return os.str();
}

/*
string cone(float radius, float height, int slices, int stacks) {
    return frustum(radius, 0.0f, height, slices, stacks);
}

string cylinder(float radius, float height, int slices, int stacks) {
    return frustum(radius, radius, height, slices, stacks);
}
*/

string torus(float innerRadius, float outerRadius, int nsides, int nrings) {
    int i, j, nPoints = 0;
    float deltaPhi, deltaTheta;
    float phi, theta, nextPhi, nextTheta;
    float cosPhi, cosTheta, sinPhi, sinTheta;
    float cosNextPhi, cosNextTheta, sinNextPhi, sinNextTheta;
    ostringstream os;

    deltaPhi = 2.0f * M_PI / nrings;
    deltaTheta = 2.0f * M_PI / nsides;
    for (i = 0; i < nrings; i++) {
        phi = i * deltaPhi;
        nextPhi = (i + 1) * deltaPhi;
        cosPhi = cosf(phi);
        sinPhi = sinf(phi);
        cosNextPhi = cosf(nextPhi);
        sinNextPhi = sinf(nextPhi);

        for (j = 0; j < nsides; j++, nPoints += 6) {
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

string torusNormals(float innerRadius, float outerRadius, int nsides, int nrings) {
    int i, j, nPoints = 0;
    float deltaPhi, deltaTheta;
    float phi, theta, nextPhi, nextTheta;
    float cosPhi, cosTheta, sinPhi, sinTheta;
    float cosNextPhi, cosNextTheta, sinNextPhi, sinNextTheta;
    ostringstream os;

    deltaPhi = 2.0f * M_PI / nrings;
    deltaTheta = 2.0f * M_PI / nsides;
    for (i = 0; i < nrings; i++) {
        phi = i * deltaPhi;
        nextPhi = (i + 1) * deltaPhi;
        cosPhi = cosf(phi);
        sinPhi = sinf(phi);
        cosNextPhi = cosf(nextPhi);
        sinNextPhi = sinf(nextPhi);

        for (j = 0; j < nsides; j++, nPoints += 6) {
            //float dXZ, nextDXZ;

            theta = j * deltaTheta;
            nextTheta = (j + 1) * deltaTheta;
            cosTheta = cosf(theta);
            sinTheta = sinf(theta);
            cosNextTheta = cosf(nextTheta);
            sinNextTheta = sinf(nextTheta);
            //dXZ = outerRadius + innerRadius * cosTheta;
            //nextDXZ = outerRadius + innerRadius * cosNextTheta;

            os << (cosTheta * sinPhi) << ' ' << sinTheta << ' ' << (cosTheta * cosPhi) << '\n';
            os << (cosNextTheta * sinNextPhi) << ' ' << sinNextTheta << ' ' << (cosNextTheta * cosNextPhi) << '\n';
            os << (cosNextTheta * sinPhi) << ' ' << sinNextTheta << ' ' << (cosNextTheta * cosPhi) << '\n';

            os << (cosTheta * sinPhi) << ' ' << sinTheta << ' ' << (cosTheta * cosPhi) << '\n';
            os << (cosTheta * sinNextPhi) << ' ' << sinTheta << ' ' << (cosTheta * cosNextPhi) << '\n';
            os << (cosNextTheta * sinNextPhi) << ' ' << sinNextTheta << ' ' << (cosNextTheta * cosNextPhi) << '\n';
        }
    }
    return os.str();
}

string torusTexCoords(int nsides, int nrings) {
    int i, j;
    ostringstream os;
    float deltaS = 1.0f / nrings;
    float deltaT = 1.0f / nsides;

    for (i = 0; i < nrings; i++) {
        float s = i * deltaS;

        for (j = 0; j < nsides; j++) {
            float t = j * deltaT;

            os << s << ' ' << t << '\n';
            os << (s + deltaS) << ' ' << (t + deltaT) << '\n';
            os << s << ' ' << (t + deltaT) << '\n';

            os << s << ' ' << t << '\n';
            os << (s + deltaS) << ' ' << t << '\n';
            os << (s + deltaS) << ' ' << (t + deltaT) << '\n';
        }
    }
    return os.str();
}

int annulusGenerator(int argc, char* argv[]) {
    float dist, smj, smn;
    int slices;

    ofstream outfile;
    dist = atof(argv[0]);
    smj = atof(argv[1]);
    smn = atof(argv[2]);
    slices = atoi(argv[3]);
    if (dist <= 0.0f || smj <= 0.0f || smn <= 0.0f || slices <= 0) {
        fputs("Error: All parameters of the annulus must be positive numbers\n", stderr);
        return 1;
    }

    outfile.open(argv[4]);
    if (!outfile.is_open()) {
        perror("ofstream.open");
        return 1;
    }
    outfile << annulus(dist, smj, smn, slices);
    outfile << annulusNormals(slices);
    outfile << annulusTexCoords(slices);
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

    if (xDim <= 0.0f || yDim <= 0.0f || zDim <= 0.0f || divisions <= 0) {
        fputs("Error: All parameters of the box must be positive numbers\n", stderr);
        return 1;
    }

    outfile.open((argc == 5) ? argv[4] : argv[3]);
    if (!outfile.is_open()) {
        perror("ofstream.open");
        return 1;
    }
    outfile << box(xDim, yDim, zDim, divisions);
    outfile << boxNormals(divisions);
    outfile << boxTexCoords(divisions);
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
    if (radius <= 0.0f || height <= 0.0f || slices <= 0 || stacks <= 0) {
        fputs("Error: All parameters of the cone must be positive numbers\n", stderr);
        return 1;
    }
    outfile.open(argv[4]);
    if (!outfile.is_open()) {
        perror("ofstream.open");
        return 1;
    }
    outfile << frustum(radius, 0.0f, height, slices, stacks);
    outfile << frustumNormals(radius, 0.0f, height, slices, stacks);
    outfile << frustumTexCoords(radius, 0.0f, slices, stacks);
    outfile.close();
    return 0;
}

int cylinderGenerator(int argc, char *argv[]) {
    float radius, height;
    int slices, stacks;
    ofstream outfile;

    radius = atof(argv[0]);
    height = atof(argv[1]);
    slices = atoi(argv[2]);
    stacks = atoi(argv[3]);
    if (radius <= 0.0f || height <= 0.0f || slices <= 0 || stacks <= 0) {
        fputs("Error: All parameters of the cylinder must be positive numbers\n", stderr);
        return 1;
    }
    outfile.open(argv[4]);
    if (!outfile.is_open()) {
        perror("ofstream.open");
        return 1;
    }
    outfile << frustum(radius, radius, height, slices, stacks);
    outfile << frustumNormals(radius, radius, height, slices, stacks);
    outfile << frustumTexCoords(radius, radius, slices, stacks);
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
    if (xRadius <= 0.0f || yRadius <= 0.0f || zRadius <= 0.0f || slices <= 0 || stacks <= 0) {
        fputs("Error: All parameters of the cylinder must be positive numbers\n", stderr); // consider creating a macro for this kind of error message
        return 1;
    }
    outfile.open(argv[5]);
    if (!outfile.is_open()) {
        perror("ofstream.open");
        return 1;
    }
    outfile << ellipsoid(xRadius, yRadius, zRadius, slices, stacks);
    outfile << ellipsoidNormals(xRadius, yRadius, zRadius, slices, stacks);
    outfile << ellipsoidTexCoords(slices, stacks);
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
    if (*endptr1 != '\0' || *endptr2 != '\0' || baseRadius < 0.0f || topRadius < 0.0f) {
        fputs("Error: The base radius and the top radius of the frustum must be non-negative numbers\n", stderr);
        return 1;
    }
    if (height <= 0 || slices <= 0 || stacks <= 0) {
        fputs("Error: Height, slices and stacks must all be positive\n", stderr);
        return 1;
    }

    outfile.open(argv[5]);
    if (!outfile.is_open()) {
        perror("ofstream.open");
        return 1;
    }
    outfile << frustum(baseRadius, topRadius, height, slices, stacks);
    outfile << frustumNormals(baseRadius, topRadius, height, slices, stacks);
    outfile << frustumTexCoords(baseRadius, topRadius, slices, stacks);
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
    if (xDim <= 0.0f || zDim <= 0.0f || divisions <= 0) {
        fputs("Error: All parameters of the plane must be positive numbers\n", stderr);
        return 1;
    }
    (argc == 4) ? outfile.open(argv[3]) : outfile.open(argv[2]);

    if (!outfile.is_open()) {
        perror("ofstream.open");
        return 1;
    }
    outfile << planeXZ(xDim, 0, zDim, divisions, UPWARDS);
    outfile << planeXZnormals(divisions, UPWARDS);
    outfile << planeTexCoords(divisions, UPWARDS);
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
    if (radius <= 0.0f || slices <= 0 || stacks <= 0) {
        fputs("Error: All parameters of the sphere must be positive numbers\n", stderr);
        return 1;
    }

    outfile.open(argv[3]);
    if (!outfile.is_open()) {
        perror("ofstream.open");
        return 1;
    }
    outfile << ellipsoid(radius, radius, radius, slices, stacks);
    outfile << ellipsoidNormals(radius, radius, radius, slices, stacks);
    outfile << ellipsoidTexCoords(slices, stacks);
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
    if (innerRadius <= 0.0f || outerRadius <= 0.0f || sides <= 0 || rings <= 0) {
        fputs("Error: All parameters of the torus must be positive numbers\n", stderr);
        return 1;
    }
    outfile.open(argv[4]);
    if (!outfile.is_open()) {
        perror("ofstream.open");
        return 1;
    }
    outfile << torus(innerRadius, outerRadius, sides, rings);
    outfile << torusNormals(innerRadius, outerRadius, sides, rings);
    outfile << torusTexCoords(sides, rings);
    outfile.close();
    return 0;
}

void getBezierPatchPoint(float u, float v, point *pv, float *res, float *pNormal, float *texCoords) {
    float dU[3];
    float dV[3];
    /* Setup */
    float m[4][4] = {
        { -1.0f,  3.0f, -3.0f, 1.0f },
        {  3.0f, -6.0f,  3.0f, 0.0f },
        { -3.0f,  3.0f,  0.0f, 0.0f },
        {  1.0f,  0.0f,  0.0f, 0.0f }
    };

    /* transpose of m equals m */
    float Px[4][4] = {
        { pv[0].x, pv[1].x, pv[2].x, pv[3].x },
        { pv[4].x, pv[5].x, pv[6].x, pv[7].x },
        { pv[8].x, pv[9].x, pv[10].x, pv[11].x },
        { pv[12].x, pv[13].x, pv[14].x, pv[15].x }
    };

    float Py[4][4] = {
        { pv[0].y, pv[1].y, pv[2].y, pv[3].y },
        { pv[4].y, pv[5].y, pv[6].y, pv[7].y },
        { pv[8].y, pv[9].y, pv[10].y, pv[11].y },
        { pv[12].y, pv[13].y, pv[14].y, pv[15].y }
    };

    float Pz[4][4] = {
        { pv[0].z, pv[1].z, pv[2].z, pv[3].z },
        { pv[4].z, pv[5].z, pv[6].z, pv[7].z },
        { pv[8].z, pv[9].z, pv[10].z, pv[11].z },
        { pv[12].z, pv[13].z, pv[14].z, pv[15].z }
    };

    float U[4] = {u * u * u, u * u, u, 1};
    float UD[4] = {3 * u * u, 2 * u, 1, 0};
    float V[4] = {v * v * v, v * v, v, 1};
    float VD[4] = {3 * v * v, 2 * v, 1, 0};

    float MdV[4];
    float MV[4];
    multMatrixVector((float *)m, V, MV);
    multMatrixVector((float *)m, VD, MdV);

    float dUM[4];
    float UM[4];
    multVectorMatrix(U, (float *)m, UM);
    multVectorMatrix(UD, (float *)m, dUM);

    float UMP[3][4];
    multVectorMatrix(UM, (float *)Px, UMP[0]);
    multVectorMatrix(UM, (float *)Py, UMP[1]);
    multVectorMatrix(UM, (float *)Pz, UMP[2]);

    float dUMP[3][4];
    multVectorMatrix(dUM, (float *)Px, dUMP[0]);
    multVectorMatrix(dUM, (float *)Py, dUMP[1]);
    multVectorMatrix(dUM, (float *)Pz, dUMP[2]);

    for (int j = 0; j < 3; j++) {
        res[j] = 0.0f;
        dU[j] = 0.0f;
        dV[j] = 0.0f;

        for (int i = 0; i < 4; i++) {
            res[j] += MV[i] * UMP[j][i];
            dU[j] += MV[i] * dUMP[j][i];
            dV[j] += MdV[i] * UMP[j][i];
        }
    }
    normalize(dU);
    normalize(dV);
    cross(dV, dU, pNormal);
    normalize(pNormal);

    texCoords[0] = u;
    texCoords[1] = v;
}

int bezierPatchGenerator(char *outfile, int tesselationLevel) {
    int p = 0;
    point pv[16];
    int divs = tesselationLevel; // change this to change the tesselation level
    ostringstream pointStr, normalsStr, texCoordsStr;
    ofstream out;

    out.open(outfile);
    if (!out.is_open()) {
        perror("ofstream.open");
        return 1;
    }
    for (int i = 0; i < patches; i++) {
        for (int j = 0; j < 16; j++) {
            pv[j] = cpoints[indexes[i][j]];
        }
        for (int u = 0; u < divs; u++) {
            float resP1[3], pNormal1[3], texCoords1[2];
            float resP2[3], pNormal2[3], texCoords2[2];
            float resP3[3], pNormal3[3], texCoords3[2];
            float resP4[3], pNormal4[3], texCoords4[2];

            for (int v = 0; v < divs; v++) {
                getBezierPatchPoint(u / (float)divs, v / (float)divs, pv, resP1, pNormal1, texCoords1);
                getBezierPatchPoint((u + 1) / (float)divs, v / (float)divs, pv, resP2, pNormal2, texCoords2);
                getBezierPatchPoint(u / (float)divs, (v + 1) / (float)divs, pv, resP3, pNormal3, texCoords3);
                getBezierPatchPoint((u + 1) / (float)divs, (v + 1) / (float)divs, pv, resP4, pNormal4, texCoords4);

                pointStr << resP1[0] << ' ' << resP1[1] << ' ' << resP1[2] << '\n';
                pointStr << resP3[0] << ' ' << resP3[1] << ' ' << resP3[2] << '\n';
                pointStr << resP4[0] << ' ' << resP4[1] << ' ' << resP4[2] << '\n';

                pointStr << resP2[0] << ' ' << resP2[1] << ' ' << resP2[2] << '\n';
                pointStr << resP1[0] << ' ' << resP1[1] << ' ' << resP1[2] << '\n';
                pointStr << resP4[0] << ' ' << resP4[1] << ' ' << resP4[2] << '\n';

                normalsStr << pNormal1[0] << ' ' << pNormal1[1] << ' ' << pNormal1[2] << '\n';
                normalsStr << pNormal3[0] << ' ' << pNormal3[1] << ' ' << pNormal3[2] << '\n';
                normalsStr << pNormal4[0] << ' ' << pNormal4[1] << ' ' << pNormal4[2] << '\n';

                normalsStr << pNormal2[0] << ' ' << pNormal2[1] << ' ' << pNormal2[2] << '\n';
                normalsStr << pNormal1[0] << ' ' << pNormal1[1] << ' ' << pNormal1[2] << '\n';
                normalsStr << pNormal4[0] << ' ' << pNormal4[1] << ' ' << pNormal4[2] << '\n';

                texCoordsStr << texCoords1[0] << ' ' << texCoords1[1] << '\n';
                texCoordsStr << texCoords3[0] << ' ' << texCoords3[1] << '\n';
                texCoordsStr << texCoords4[0] << ' ' << texCoords4[1] << '\n';

                texCoordsStr << texCoords2[0] << ' ' << texCoords2[1] << '\n';
                texCoordsStr << texCoords1[0] << ' ' << texCoords1[1] << '\n';
                texCoordsStr << texCoords4[0] << ' ' << texCoords4[1] << '\n';
                
                p += 6;
            }
        }
    }
    out << (to_string(p) + "\n" + pointStr.str() + normalsStr.str() + texCoordsStr.str());
    out.close();
    
    return 0;
}

int bezierPatchParser(char *patch) {
    int lIndex;
    int i, j;
    char line[BUFF_SIZE];
    FILE *f = fopen(patch, "r");

    if (!f)
        return 1;

    fscanf(f, "%d\n", &lIndex);
    patches = lIndex;
    printf("%d\n", patches);
    indexes = (int **) malloc(sizeof(int *) * lIndex);
    if (!indexes)
        return 1;

    for (i = 0; i < lIndex; i++) {
        indexes[i] = (int *) malloc(sizeof(int) * 16);
        if (!indexes[i]) {
            // free previously allocated memory before returning non zero
            for (--i; i >= 0; --i) {
                free(indexes[i]);
            }
            free(indexes);

            return 1;
        }
        memset(line, 0, BUFF_SIZE);
        fgets(line, BUFF_SIZE, f);
        char* ind = NULL;
        for(j = 0, ind = strtok(line,", "); ind && j < 16; ind = strtok(NULL, ", "), j++)
            indexes[i][j] = atoi(ind);
    }
    fscanf(f, "%d\n", &ncpoints);

    cpoints = (point *) malloc(sizeof(point) * ncpoints);
    if (!cpoints) {
        for (i = 0; i < lIndex; i++) {
            free(indexes[i]);
        }
        free(indexes);

        return 1;
    }

    for (i = 0; i < ncpoints; i++) {
        memset(line, 0, BUFF_SIZE);
        fgets(line, BUFF_SIZE, f);
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
          "| patch       | tesselationLevel patchFile                |\n"
          "+-------------+-------------------------------------------+\n"
          , stream
         );
}

int main(int argc, char *argv[]) {
    int rval;

    if (argc == 1) {
        usage(argv[0], stderr);
        return 1;
    }

    string primitive(argv[1]);
    if (primitive == "annulus" && argc == 7) {
        rval = annulusGenerator(argc - 2, &argv[2]);
    }
    else if (primitive == "box" && (argc == 6 || argc == 7)) {
        rval = boxGenerator(argc - 2, &argv[2]);
    }
    else if (primitive == "cone" && argc == 7) {
        rval = coneGenerator(argc - 2, &argv[2]);
    }
    else if (primitive == "cylinder" && argc == 7) {
        rval = cylinderGenerator(argc - 2, &argv[2]);
    }
    else if (primitive == "ellipsoid" && argc == 8) {
        rval = ellipsoidGenerator(argc - 2, &argv[2]);
    }
    else if (primitive == "frustum" && argc == 8) {
        rval = frustumGenerator(argc - 2, &argv[2]);
    }
    else if (primitive == "plane" && (argc == 5 || argc == 6)) {
        rval = planeGenerator(argc - 2, &argv[2]);
    }
    else if (primitive == "sphere" && argc == 6) {
        rval = sphereGenerator(argc - 2, &argv[2]);
    }
    else if (primitive == "torus" && argc == 7) {
        rval = torusGenerator(argc - 2, &argv[2]);
    }
    else if (primitive == "patch" && argc == 5) {
        rval = bezierPatchParser(argv[3]);
        if (!rval)
            rval = bezierPatchGenerator(argv[4], atoi(argv[2]));
    }
    else if (primitive == "--help") {
        usage(argv[0], stdout);
        rval = 0;
    } else {
        usage(argv[0], stderr);
        rval = 1;
    }
    return rval;
}
