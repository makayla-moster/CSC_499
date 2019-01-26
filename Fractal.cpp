//g++ -w -o Makayla.exe gl_utils.cpp maths_funcs.cpp Fractal.cpp libglew32.dll.a libglfw3dll.a -I include -lglfw3 -lgdi32 -lopengl32 -L ./ -lglew32 -lglfw3 //Not using
//g++ -w -o Makayla.exe Fractal.cpp libglew32.dll.a libglfw3dll.a -I include -lOpenGL32 -L ./ -lglew32 -lglfw3 																							// Not using
//g++ -w -o Makayla.exe gl_utils.cpp maths_funcs.cpp Fractal.cpp libglfw3dll.a libglew32.dll.a -I include -lglfw3 -lgdi32 -lopengl32												// Not using
// g++ -w -o Makayla.exe gl_utils.cpp maths_funcs.cpp Fractal.cpp libglfw3dll.a libglew32.dll.a -I include -lgdi32 -lopengl32 -L ./ -lglew32 -lglfw3				// Not using

// g++ -w -o monkm.exe gl_utils.cpp maths_funcs.cpp Fractal.cpp libglfw3dll.a libglew32.dll.a -I include -lgdi32 -lopengl32 -L ./ -lglew32 -lglfw3

#include "gl_utils.h"
#include "maths_funcs.h"
#include <GL/glew.h>																														/* include GLEW and new version of GL on Windows */
#include <GLFW/glfw3.h>     																										/* GLFW helper library */
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <stack>
#include <vector>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <time.h>
#include <assert.h>
#include <string.h>
#include <stdarg.h>
#define GL_LOG_FILE "gl.log"
using namespace std;

void multiplyNew(GLfloat matrix1[], GLfloat matrix2[], GLfloat result[]){  			// multiply for the user interaction

	result[0] = (matrix1[0]*matrix2[0])+(matrix1[4]*matrix2[1])+(matrix1[8]*matrix2[2])+(matrix1[12]*matrix2[3]);
	result[4] = (matrix1[0]*matrix2[4])+(matrix1[4]*matrix2[5])+(matrix1[8]*matrix2[6])+(matrix1[12]*matrix2[7]);
	result[8] = (matrix1[0]*matrix2[8])+(matrix1[4]*matrix2[9])+(matrix1[8]*matrix2[10])+(matrix1[12]*matrix2[11]);
	result[12] = (matrix1[0]*matrix2[12])+(matrix1[4]*matrix2[13])+(matrix1[8]*matrix2[14])+(matrix1[12]*matrix2[15]);
	result[1] = (matrix1[1]*matrix2[0])+(matrix1[5]*matrix2[1])+(matrix1[9]*matrix2[2])+(matrix1[13]*matrix2[3]);
	result[5] = (matrix1[1]*matrix2[4])+(matrix1[5]*matrix2[5])+(matrix1[9]*matrix2[6])+(matrix1[13]*matrix2[7]);
	result[9] = (matrix1[1]*matrix2[8])+(matrix1[5]*matrix2[9])+(matrix1[9]*matrix2[10])+(matrix1[13]*matrix2[11]);
	result[13] = (matrix1[1]*matrix2[12])+(matrix1[5]*matrix2[13])+(matrix1[9]*matrix2[14])+(matrix1[13]*matrix2[15]);
	result[2] = (matrix1[2]*matrix2[0])+(matrix1[6]*matrix2[1])+(matrix1[10]*matrix2[2])+(matrix1[14]*matrix2[3]);
	result[6] = (matrix1[2]*matrix2[4])+(matrix1[6]*matrix2[5])+(matrix1[10]*matrix2[6])+(matrix1[14]*matrix2[7]);
	result[10] = (matrix1[2]*matrix2[8])+(matrix1[6]*matrix2[9])+(matrix1[10]*matrix2[10])+(matrix1[14]*matrix2[11]);
	result[14] = (matrix1[2]*matrix2[12])+(matrix1[6]*matrix2[13])+(matrix1[10]*matrix2[14])+(matrix1[14]*matrix2[15]);
	result[3] = (matrix1[3]*matrix2[0])+(matrix1[7]*matrix2[1])+(matrix1[11]*matrix2[2])+(matrix1[15]*matrix2[3]);
	result[7] = (matrix1[3]*matrix2[4])+(matrix1[7]*matrix2[5])+(matrix1[11]*matrix2[6])+(matrix1[15]*matrix2[7]);
	result[11] = (matrix1[3]*matrix2[8])+(matrix1[7]*matrix2[9])+(matrix1[11]*matrix2[10])+(matrix1[15]*matrix2[11]);
	result[15] = (matrix1[3]*matrix2[12])+(matrix1[7]*matrix2[13])+(matrix1[11]*matrix2[14])+(matrix1[15]*matrix2[15]);
}

void multiply(GLfloat matrix1[], GLfloat matrix2[], GLfloat result1[]){   			// this is to multiply a 4x4 and a 4x1 matrix together
	result1[0] = (matrix1[0]*matrix2[0]) + (matrix1[4]*matrix2[1]) + (matrix1[8]*matrix2[2]) + (matrix1[12]*matrix2[3]);
	result1[1] = (matrix1[1]*matrix2[0]) + (matrix1[5]*matrix2[1]) + (matrix1[9]*matrix2[2]) + (matrix1[13]*matrix2[3]);
	result1[2] = (matrix1[2]*matrix2[0]) + (matrix1[6]*matrix2[1]) + (matrix1[10]*matrix2[2]) + (matrix1[14]*matrix2[3]);
	result1[3] = (matrix1[3]*matrix2[0]) + (matrix1[7]*matrix2[1]) + (matrix1[11]*matrix2[2]) + (matrix1[15]*matrix2[3]);
}

GLfloat* multiplyAgain(GLfloat matrix1[], GLfloat matrix2[], GLfloat result1[]){
	for (int i = 0; i < 16; i++){
		for (int j = 0; j < 4; j++){
			result1[i] = result1[i] + (matrix2[(i / 4) * 4 + j] * matrix1[(i % 4)+(4*j)]);
		}
	}
	return result1;
}

string generatePattern(){																												//Generates a pattern to create a tree.
    int numIts = 4; 																														// Number of iterations
    string pattern = "F"; //"[X]";    																					// Using F for the pattern

    for (int i = 0; i < numIts; i++){
        string newPattern = "";
        for (int idx = 0; idx < pattern.length(); idx++){
            if (pattern.substr(idx,1).compare("F") == 0)
			newPattern += "F[F][-F][+F]";   //"F[-F][F][+F]"
            else if (pattern.substr(idx,1).compare("X") == 0)
                newPattern += "F-[[X]+X]+F[+FX]-X";
            else{
                newPattern += pattern.substr(idx,1);
            }
        }
        pattern = newPattern;
    }

    return pattern;
}

int countF(string pattern){
	int count = 1;
	for (int idx = 0; idx < pattern.length(); idx++){															// Loops through string pattern to find number of 'F' within the pattern.
		if (pattern.substr(idx, 1).compare("F") == 0){
			count ++;
		}
	}
	return count;
}


int countbracket(string pattern) {
	int countBracket = 0;
	for (int idx = 0; idx < pattern.length(); idx++){															// Loops through string pattern to find number of ']' within the pattern.
		if (pattern.substr(idx, 1).compare("]") == 0){
			countBracket ++;
		}
	}

	return countBracket;
}

int countLbracket(string pattern) {
	int countLBracket = 0;
	for (int idx = 0; idx < pattern.length(); idx++){															// Loops through string pattern to find number of ']' within the pattern.
		if (pattern.substr(idx, 1).compare("[]") == 0){
			countLBracket ++;
		}
	}

	return countLBracket;
}

int countLabel(string modelName, char label[]){																	// Counts number of labels within an OBJ file
	int numLab = 0;

	FILE *objFile;
	objFile = fopen(modelName.c_str(),"r");

	char buf[128];
	while (fscanf(objFile, "%s", &buf) != EOF){
	    if (strcmp(buf,label) == 0)
		numLab++;
	}

	cout << "Model has " << numLab << " " << label << "\n";
	fclose(objFile);
	return numLab;
}

void loadVertices(string modelName, GLfloat verts[]){														// Loads the vertices of an OBJ file into the program
	cout << "Loading vertices\n";
	int numVert = 0;

	FILE *objFile;
	objFile = fopen(modelName.c_str(),"r");

	float maxX = -100000000.0;
	float minX = 100000000.0;
	float maxY = -100000000.0;
	float minY = 100000000.0;
	float maxZ = -100000000.0;
	float minZ = 100000000.0;

	char buf[128];
	char label[] = "v";
	float a, b, c;
	while (fscanf(objFile, "%s", &buf) != EOF){
	    if (strcmp(buf,label) == 0){
		fscanf(objFile, "%f %f %f\n", &a, &b, &c);
		if (a > maxX)
		    maxX = a;
		if (a < minX)
		    minX = a;
		if (b > maxY)
		    maxY = b;
		if (b < minY)
		    minY = b;
		if (c > maxZ)
		    maxZ = c;
		if (c < minZ)
		    minZ = c;
		verts[3*numVert + 0] = 1.0*a;
		verts[3*numVert + 1] = 1.0*b;
		verts[3*numVert + 2] = 1.0*c;
		numVert++;
	    }
	}

	float scaleX = maxX-minX;
	float scaleY = maxY-minY;
	float scaleZ = maxZ-minZ;
	float transX = 0.5*(maxX+minX);
	float transY = 0.5*(maxY+minY);
	float transZ = 0.5*(maxZ+minZ);
	cout << "scales: " << scaleX << ", " << scaleY << ", " << scaleZ << endl;

	for (int i = 0; i < numVert; i++){
	    verts[3*i+0] = (verts[3*i+0] - transX)/scaleX;
	    verts[3*i+1] = (verts[3*i+1] - transY)/scaleY;
	    verts[3*i+2] = (verts[3*i+2] - transZ)/scaleZ;
	}

	fclose(objFile);
	cout << "Done loading vertices\n";
}


void computeFaceNormals(GLfloat faceNormals[], GLfloat verts[], GLint faces[], int numFaces){
	for (int i = 0; i < numFaces; i++){																						// Computes the number of face normals within an OBJ file
		int idx1 = faces[i*3 + 0];
		int idx2 = faces[i*3 + 1];
		int idx3 = faces[i*3 + 2];
		float ux = verts[idx1*3 + 0] - verts[idx2*3 + 0];
		float uy = verts[idx1*3 + 1] - verts[idx2*3 + 1];
		float uz = verts[idx1*3 + 2] - verts[idx2*3 + 2];
		float vx = verts[idx1*3 + 0] - verts[idx3*3 + 0];
		float vy = verts[idx1*3 + 1] - verts[idx3*3 + 1];
		float vz = verts[idx1*3 + 2] - verts[idx3*3 + 2];
		float nx = uy*vz - uz*vy;
		float ny = uz*vx - ux*vz;
		float nz = ux*vy - uy*vx;
		float mag = sqrt(nx*nx + ny*ny + nz*nz);
		faceNormals[3*i + 0] = nx/mag;
		faceNormals[3*i + 1] = ny/mag;
		faceNormals[3*i + 2] = nz/mag;
	}

}

void computeVertNormals(GLfloat normals[], GLfloat verts[], int numVerts, GLint faces[], int numFaces, GLfloat faceNormals[]){
	for (int i = 0; i < numVerts; i++){																						//Computes the number of vertex normals in an OBJ file
		float avgX = 0.0;
		float avgY = 0.0;
		float avgZ = 0.0;
		int numF_vert = 0;

		//find all the faces that contain this vertex
		for (int j = 0; j < numFaces; j++){
			int found = 0;
			for (int k = 0; k < 3; k++)
				if (faces[j*3 + k] == i)
					found = 1;
			if (found){
				avgX += faceNormals[j*3 + 0];
				avgY += faceNormals[j*3 + 1];
				avgZ += faceNormals[j*3 + 2];
				numF_vert++;
			}
		}

		avgX /= numF_vert;
		avgY /= numF_vert;
		avgZ /= numF_vert;
		normals[i*3 + 0] = avgX;
		normals[i*3 + 1] = avgY;
		normals[i*3 + 2] = avgZ;
	}
}

void loadFaces(string modelName, GLint faces[]){    														// Loads a Maya OBJ file into the program
    cout << "Loading new faces\n";

    FILE *objFile;
    objFile = fopen(modelName.c_str(),"r");
    int numFaces = 0;
    char buf[128];
    char label[] = "f";
    float a, b, c, a1, b1, c1, a2, b2, c2;
    while (fscanf(objFile, "%s", &buf) != EOF){
        if (strcmp(buf,label) == 0){
        fscanf(objFile, "%f/%f/%f %f/%f/%f %f/%f/%f\n", &a, &b, &c, &a1, &b1, &c1, &a2, &b2, &c2);
                faces[3*numFaces+0] = a-1;
                faces[3*numFaces+1] = a1-1;
                faces[3*numFaces+2] = a2-1;
                numFaces++;
        }
    }

    fclose(objFile);
    cout << "Done loading faces\n";
}

float RandomFloat(float a, float b) {																						// Generates a random float value between 2 float values which can be negative
    float random = ((float) rand()) / (float) RAND_MAX;
    float diff = b - a;
    float r = random * diff;
    return a + r;
}


/* Begin Code for User Interaction feature */																		// This is the code for User Interaction
																																								// This code will allow the user to rotate, skew, and translate the tree
float x, y, z, rx, ry, rz = 0;																									// Variables initialized for translation and rotational matrices
float sy = 1.0;																																	// Variable - scale matrix y component
float sx = 1.0;																																	// Variable - scale matrix x component
float sz = 0.0;																																	// Variable - scale matrix z component
float d = 1.6;																																	// Variable - initializes the skew amount
float fov = 67 * 3.14159 /180.0;																								// Variable - initilizes field of view
float aspect = 1.0;																															// Variable - initializes aspect ratio
float near = 0.01;																															// Variable - initializes near
float far = 100.0;																															// Variable - initializes far
float range = tan(fov*0.5)*near;																								// Variable - initializes the range of view
float Sx = (2*near)/((range*aspect)+(range*aspect));														// Variable - for proMat
float Sy = near/range;																													// Variable - for proMat
float Sz = -(far+near)/(far-near);																							// Variable - for proMat
float Pz = -(2*far*near)/(far-near);																						// Variable - for proMat

GLfloat proMat[] = {Sx, 0.0f, 0.0f, 0.0f,																				// Matrix - ??
				0.0f, Sy, 0.0f, 0.0f,
				0.0f, 0.0f, Sz, -1.0,
				0.0f, 0.0f, Pz, 1.0f};

GLfloat ortho[] =  {1,0,0,0,																										// Matrix - Orthogonal viewing matrix
							0,1,0,0,
							0,0,0,0,
							0,0,0,1};

GLfloat lookAt[] = {1.0f, 0.0f, -0.0f, 0.0f,																		// Matrix - Look At matrix
					0.0f, 1.0f, -0.0f, 0.0f,
					0.0f, 0.0f, 1.0f, 0.0f,
					-0.0f, -0.0f, -1.0f, 1.0f};

GLfloat view[] =  { 1,0,0,0,																										// Matrix - View matrix
					0,1,0,0,
					0,0,1,0,
					0,0,0,1};
GLfloat* viewResult = new float[16];																						// Matrix - created for holding two multiplied matrices
GLfloat* viewMat = new float[16];																								// Matrix - created for holding two multiplied matrices

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods){ // User interaction portion
	if(key == GLFW_KEY_A && action == GLFW_PRESS){
		cout<<"A pressed\n"; 																												// Moves the tree left
		x -= 0.1;
	}
	if(key == GLFW_KEY_D && action == GLFW_PRESS){
		cout<<"D pressed\n"; 																												// Moves the tree right
		x += 0.1;
	}
	if(key == GLFW_KEY_W && action == GLFW_PRESS){
		cout<<"W pressed\n"; 																												// Moves the tree upward
		y += 0.1;
	}
	if(key == GLFW_KEY_S && action == GLFW_PRESS){
		cout<<"S pressed\n"; 																												// Moves the tree downward
		y -= 0.1;
	}
	if(key == GLFW_KEY_X && action == GLFW_PRESS){
		cout<<"X pressed\n";																												// Moves the tree positively along the Z axis
		z += 0.1;
	}
	if(key == GLFW_KEY_Z && action == GLFW_PRESS){
		cout<<"Z pressed\n";																												// Moves the tree negatively along the Z axis
		z -= 0.1;
	}
	// if(key == GLFW_KEY_I && action == GLFW_PRESS){
	// 	cout<<"I pressed\n";																												// Scales the tree's x value negatively
	// 	sx -= 0.1;
	// 	//cout<<"-sx :"<<sx<<"\n";
	// }
	// if(key == GLFW_KEY_K && action == GLFW_PRESS){
	// 	cout<<"K pressed\n";																												// Scales the tree's x value positively
	// 	sx += 0.1;
	// 	//cout<<"+sx :"<<sx<<"\n";
	// }
	// if(key == GLFW_KEY_L && action == GLFW_PRESS){
	// 	cout<<"L pressed\n";																												// Scales the tree's y value positively
	// 	sy += 0.1;
	// 	//cout<<"+sy :"<<sy<<"\n";
	// }
	// if(key == GLFW_KEY_J && action == GLFW_PRESS){
	// 	cout<<"J pressed\n";																												// Scales the tree's y value negatively
	// 	sy -= 0.1;
	// 	//cout<<"-sy :"<<sy<<"\n";
	// }
	// if(key == GLFW_KEY_M && action == GLFW_PRESS){
	// 	cout<<"M pressed\n";																												// Scales the tree's z value positively
	// 	sz += 0.1;
	// 	//cout<<"+sz :"<<sz<<"\n";
	// }
	// if(key == GLFW_KEY_N && action == GLFW_PRESS){
	// 	cout<<"N pressed\n";																												// Scales the tree's z value negatively
	// 	sz -= 0.1;
	// 	//cout<<"-sz :"<<sz<<"\n";
	// }
	if(key == GLFW_KEY_R && action == GLFW_PRESS){
		cout<<"R pressed\n";																												// Rotates the tree negatively along the X axis
		rx -= 0.1;
	}
	if(key == GLFW_KEY_T && action == GLFW_PRESS){
		cout<<"T pressed\n";																												// Rotates the tree positively along the X axis
		rx += 0.1;
	}
	if(key == GLFW_KEY_F && action == GLFW_REPEAT){
		cout<<"F pressed\n";																												// Rotates the tree negatively along the Y axis
		ry -= 0.1;
	}
	if(key == GLFW_KEY_G && action == GLFW_REPEAT){
		cout<<"G pressed\n";																												// Rotates the tree positively along the Y axis
		ry += 0.1;
	}
	if(key == GLFW_KEY_V && action == GLFW_PRESS){
		cout<<"V pressed\n";																												// Rotates the tree negatively along the Z axis
		rz -= 0.1;
		//cout<<"rz + :"<<rz<<"\n";
	}
	if(key == GLFW_KEY_B && action == GLFW_PRESS){
		cout<<"B pressed\n";																												// Rotates the tree positively along the Z axis
		rz += 0.1;
		//cout<<"rz - :"<<rz<<"\n";
	}
	if(key == GLFW_KEY_Y && action == GLFW_PRESS){
		cout<<"Y pressed\n";																												// Changes the field of view positively
		fov += 0.1;
		//cout<<d<<"\n";
	}
	if(key == GLFW_KEY_U && action == GLFW_PRESS){
		cout<<"U pressed\n";																												// Changes the field of view negatively
		fov -= 0.1;
		//cout<<d<<"\n";
	}
	if(key == GLFW_KEY_E && action == GLFW_PRESS){
		cout<<"E pressed\n";																												// Changes the aspect positively ??
		aspect += 0.1;
		//cout<<d<<"\n";
	}
	if(key == GLFW_KEY_C && action == GLFW_PRESS){
		cout<<"C pressed\n";																												// Changes the aspect negatively ??
		aspect -= 0.1;
		//cout<<d<<"\n";
	}
	if(key == GLFW_KEY_H && action == GLFW_PRESS){
		cout << "H pressed - changed to promat\n";																	// Changing the viewing matrix to the proMat
		viewMat = proMat;
	}
}		//Adds in user interaction feature

/* End code for user interaction feature */

int g_gl_width = 640;																														// Variable - width of the window displayed
int g_gl_height = 480;																													// Variable - height of the window displayed
GLFWwindow* g_window = NULL;

int main() {																																		// MAIN FUNCTION WHERE CODE WILL BE RUN
	const GLubyte *renderer;
	const GLubyte *version;
	GLuint vao;																																		// Variable - Vertex Array Object initialized
	GLuint vbo;																																		// Variable - Vertex Buffer Object initialized

	int rotation;																																	// Variable - Rotation of branches.
	int count;																																		// Variable - Counts number of 'F' in string.
	int countBracket;																															// Variable - Counts number of ']' in string.
	int countLBracket;

	string pattern = generatePattern();																						// Variable - Generates string pattern to make tree from.
	// cout << pattern << endl << endl;
	count = countF(pattern);																											// Variable - Function to count the number of 'F' in the string.
	countBracket = countbracket(pattern);																					// Variable - Function to count the number of ']' in the string.
	int totalCount = count + countBracket;																				// Variable - Total amount of points, including the backtracking points that are added for the lines.
	int totalLeafCount = count*3 + countBracket;																	// Variable - Total amount of points*3 to create leaves
	// cout << totalCount << endl;

	GLfloat branchPoints[totalCount*3]; 																					//List of points to make the branches. Includes extra points for lines.
	GLfloat leafPoints[totalLeafCount*3];																					//List of points to place the leaves - only on the ends of branches though.
	GLfloat* result1 = new float[4];																							//This is a new 4x1 matrix to acquire the new heading.
	stack<float> PositionStack;																										//Stack to put branch point positions on.
	stack<float> HeadingStack;																										//Stack to put branch headings on.
	int pointsCount = 0;																													//Counts number of points needed to make a matrix of points.
	int leafCount = 0;																														//Counts number of points needed to make a matrix for leaves.

	float z3;
	float numRotateZ = (90 * 3.14159) / 180;																			// Variable - For the first rotation, so the trunk is 90 degrees from the bottom of the screen.
	float numRotateZ2 = 0; //(30 * 3.14159) / 180;																// Variable - For rotation of leaf.
	float leafRotation = 0;																												// Variable - ??
	float numRotateX = -(90 * 3.14159) / 180;																			// Variable - Rotates the leaf OBJ to be upright in radians.
	float numRotateY = (0 * 3.14159) / 180;																				// Variable - Rotates the leaf along the y-axis.
	GLfloat dx = 0;																																// Variable - For leaf translation along the x-axis.
	GLfloat dy = 0;																																// Variable - For leaf translation along the y-axis.
	GLfloat dz = 0;																																// Variable - For leaf translation along the z-axis.
	GLfloat scaleXNum = 1; //.0095;																								// Variable - Scales the leaf.
	GLfloat scaleYNum = 1; //.1095;																								// Variable - Scales the leaf.
	GLfloat scaleZNum = 1; //.1095;																								// Variable - Scales the leaf.
	GLfloat currentPosition[] = {0.0f, 0.01f, 0.0f, 1.0f};												// Vector - Beginning current position of the tree.
	GLfloat currentHeading[] = {0.0f, 0.5f, 0.0f, 0.0f};													// Vector - Beginning current heading of the tree.
	GLfloat rotateZ1[] = 																													// Matrix - Rotation matrix 1 for the z-axis.
		{cos(numRotateZ),sin(numRotateZ),0,0,
		-sin(numRotateZ),cos(numRotateZ),0,0,
		0,0,1,0,
		0,0,0,1};
	GLfloat rotateZ2[] = 																													// Matrix - Rotation matrix 2 for the z-axis.
		{cos(numRotateZ2),sin(numRotateZ2),0,0,
		-sin(numRotateZ2),cos(numRotateZ2),0,0,
		0,0,1,0,
		0,0,0,1};
	GLfloat rotateZ3[] = 																													// Matrix - Rotation matrix 3 for the z-axis.
		{cos(z3),sin(z3),0,0,
		-sin(z3),cos(z3),0,0,
		0,0,1,0,
		0,0,0,1};
	GLfloat scale1[] =																														// Matrix - Scale matrix 1.
		{scaleXNum,0,0,0,
		 0,scaleYNum,0,0,
		 0,0,scaleZNum,0,
		 0,0,0,1};
	GLfloat rotateX1[] = 																													// Matrix - Rotation matrix 1 for the x-axis.
		{1,0,0,0,
		 0,cos(numRotateX),sin(numRotateX),0,
		 0,-sin(numRotateX),cos(numRotateX),0,
		 0,0,0,1};
	GLfloat rotateY1[] =																													// Matrix - Rotation matrix 2  for the y-axis.
		{cos(numRotateY),0,-sin(numRotateY),0,
		 0,1,0,0,
		 sin(numRotateY),0,cos(numRotateY),0,
		 0,0,0,1};
	GLfloat translateMat[] =																											// Matrix - Translation matrix.
		{1,0,0,0,
		 0,1,0,0,
		 0,0,1,0,
		 dx,dy,dz,1};
	GLfloat identity[] =																													// Matrix - Identity Matrix.
		{1,0,0,0,
		 0,1,0,0,
		 0,0,1,0,
		 0,0,0,1};

/* Code to begin the User Interaction feature of this project */								// User Interaction Code - user will be able to translate, rotate, and skew the tree.

	GLfloat translate[] =																													// Matrix - translation matrix
		{1,0,0,0,
		 0,1,0,0,
		 0,0,1,0,
		 x,y,z,1};

	GLfloat scale[] =  																														// Matrix - scaling matrix
		{sx,0,0,0,
			0,sy,0,0,
			0,0,sz,0,
			0,0,0,1};

	GLfloat rotateX[] = 																													// Matrix - X Rotation Matrix
		{1,0,0,0,
		 0, cos(rx),sin(rx),0,
		 0, -sin(rx),cos(rx),0,
		 0,0,0,1};

	GLfloat rotateY[] = 																													// Matrix - Y Rotation Matrix
		{cos(ry),0,-sin(ry),0,
			0, 1,0,0,
			sin(ry),0,cos(ry),0,
			0,0,0,1};

	GLfloat rotateZ[] = 																													// Matrix - Z Rotation Matrix
		{cos(rz),sin(rz),0,0,
		-sin(rz), cos(rz),0,0,
		0, 0,1,0,
		0,0,0,1};

	GLfloat trans[] = 																														// Matrix - Translation Matrix
		{1,0,0,0,
		 0,1,0,0,
		 0,0,1,0,
	   0,0,0,1};

	GLfloat proMat[] = 																														// Matrix - ??
		{Sx, 0.0f, 0.0f, 0.0f,
		 0.0f, Sy, 0.0f, 0.0f,
	   0.0f, 0.0f, Sz, -1.0,
		 0.0f, 0.0f, Pz, 1.0f};

	GLfloat lookAt[] = 																														// Matrix - ??
		{1.0f, 0.0f, -0.0f, 0.0f,
		 0.0f, 1.0f, -0.0f, 0.0f,
		 0.0f, 0.0f, 1.0f, 0.0f,
		-0.0f, -0.0f, -1.0f, 1.0f};

	GLfloat ortho[] = 																														// Matrix - Orthogonal Viewing Matrix
		{1,0,0,0,
		 0,1,0,0,
		 0,0,0,0,
		 0,0,0,1};

	GLfloat* result = new float[16];																							// Variable - New 4x4 matrix
	GLfloat* resultRotation = new float[16];
	GLfloat* resultRotation2 = new float[16];

	int counterF = 0;
	int counter12 = 0;
	bool b1 = false;
	float box = 0.5;
	GLfloat lastPos[] = {0.0f, 0.0f, 0.0f, 0.0f};

	/* End code for user interaction feature */

	branchPoints[pointsCount + 0] = currentPosition[0];														// These lines add the first set of points to the list of points.
	branchPoints[pointsCount + 1] = currentPosition[1];
	branchPoints[pointsCount + 2] = currentPosition[2];
	pointsCount += 3;																															// Increments pointsCount by 3

	for (int idx = 0; idx < pattern.length(); idx++){															// Parser - begins by going through each character of the pattern

	// if (((currentPosition[0] <= box && currentPosition[1] <= box) && (currentPosition[1] <= box && currentPosition[2] <= box)) && ((currentPosition[0] >= -box && currentPosition[1] >= -box) && (currentPosition[1] >= -box && currentPosition[2] >= -box))){

	// if (currentPosition[0] < box && currentPosition[1] < box){

			if (pattern.substr(idx,1).compare("[") == 0){																// IF the character matches [

					PositionStack.push(currentPosition[3]);																		// Then push the currentPosition onto the PositionStack.
					PositionStack.push(currentPosition[2]);
					PositionStack.push(currentPosition[1]);
					PositionStack.push(currentPosition[0]);

					HeadingStack.push(currentHeading[3]);																			// And push the currentHeading onto the HeadingStack.
					HeadingStack.push(currentHeading[2]);
					HeadingStack.push(currentHeading[1]);
					HeadingStack.push(currentHeading[0]);

					counterF++;
					// cout << counterF << endl;
				}


			else if (pattern.substr(idx, 1).compare("]") == 0){													// ELSE IF the character matches ]
				// if (currentPosition[0] < box && currentPosition[1] < box){
					leafPoints[leafCount + 0] = currentPosition[0];														// Add the point to the leaf array (where leaves will be placed)
					leafPoints[leafCount + 1] = currentPosition[1];
					leafPoints[leafCount + 2] = currentPosition[2];
					leafCount += 3;																													// Increment leafCount by 3
				// }


					// cout << endl << "After:" << endl;
					// cout << "Position X: " << currentPosition[0] << endl;
					// cout << "Position Y: " << currentPosition[1] << endl;
					// cout << "Position Z: " << currentPosition[2] << endl << endl;
					// cout << "Heading X: " << currentHeading[0] << endl;
					// cout << "Heading Y: " << currentHeading[1] << endl;
					// cout << "Heading Z: " << currentHeading[2] << endl << endl;

					currentPosition[0] = PositionStack.top();																	// Sets the current position back to the top of the stack.
					PositionStack.pop();																											// Pops the current position from the top of the stack.
					currentPosition[1] = PositionStack.top();
					PositionStack.pop();
					currentPosition[2] = PositionStack.top();
					PositionStack.pop();
					currentPosition[3] = PositionStack.top();
					PositionStack.pop();


					branchPoints[pointsCount + 0] = currentPosition[0];												// Adds the currentPosition to the list of branching points.
					branchPoints[pointsCount + 1] = currentPosition[1];
					branchPoints[pointsCount + 2] = currentPosition[2];
					pointsCount += 3;																													// Increments pointsCount by 3

					currentHeading[0] = HeadingStack.top();																		// Sets the currentHeading to the top of the HeadingStack.
					HeadingStack.pop();																												// Pops the currentHeading from the top of the stack.
					currentHeading[1] = HeadingStack.top();
					HeadingStack.pop();
					currentHeading[2] = HeadingStack.top();
					HeadingStack.pop();
					currentHeading[3] = HeadingStack.top();
					HeadingStack.pop();

					counterF --;
				// }
			}

			else if (pattern.substr(idx, 1).compare("F") == 0){													// ELSE IF the character matches F
					// cout << endl << "Before:" << endl;
					// cout << "Position X: " << currentPosition[0] << endl;
					// cout << "Position Y: " << currentPosition[1] << endl;
					// cout << "Position Z: " << currentPosition[2] << endl << endl;
					// cout << "Heading X: " << currentHeading[0] << endl;
					// cout << "Heading Y: " << currentHeading[1] << endl;
					// cout << "Heading Z: " << currentHeading[2] << endl << endl;

					GLfloat lastPosX = currentPosition[0];																		// Save the last current Position
					GLfloat lastPosY = currentPosition[1];
					GLfloat lastPosZ = currentPosition[2];

					// if (pointsCount > 6){																											// IF pointsCount is greater than 6
					// 	float currentZ = RandomFloat(-1.0, 1.0);																// Generate a random float in between -1 and 1
					// 	currentHeading[2] = currentZ;																						// Set that float to be the z value of the current heading
					// }

					currentPosition[0] += currentHeading[0]*.15;																// Add the currentPosition and the currentHeading together
					currentPosition[1] += currentHeading[1]*.15;																// Multiply by .2 to change the height of the tree
					currentPosition[2] += currentHeading[2]*.15;
					currentPosition[3] += currentHeading[3];																	// This never changes - determines whether it is a point or a line

					// cout << "Position X: " << currentPosition[0] << endl;
					// cout << "Position Y: " << currentPosition[1] << endl;
					// cout << "Position Z: " << currentPosition[2] << endl << endl;

					GLfloat midX = (lastPosX + currentPosition[0]) / 2;												// Get the midPoint between the currentPosition and the last saved position
					GLfloat midY = (lastPosY + currentPosition[1]) / 2;
					GLfloat midZ = (lastPosZ + currentPosition[2]) / 2;

					// if (((currentPosition[0] <= box || currentPosition[1] <= box) || (currentPosition[1] <= box || currentPosition[2] <= box)) || ((currentPosition[0] >= -box || currentPosition[1] >= -box) || (currentPosition[1] >= -box || currentPosition[2] >= -box))){
					if (currentPosition[0] <= box && currentPosition[1] <= box){

						leafPoints[leafCount + 0] = midX;																					// Adds the midpoint to the leaf array
						leafPoints[leafCount + 1] = midY;
						leafPoints[leafCount + 2] = midZ;
						leafCount += 3;																														// Increments leafCount by 3

						branchPoints[pointsCount + 0] = currentPosition[0];												//Adds the currentPosition to the list of branching points
						branchPoints[pointsCount + 1] = currentPosition[1];
						branchPoints[pointsCount + 2] = currentPosition[2];
						pointsCount += 3;																													// Increments pointsCount by 3
					}
					else {
						// cout << "HERE" << endl;
						b1 = true;

						// cout << "Position X: " << currentPosition[0] << endl;
						// cout << "Position Y: " << currentPosition[1] << endl;
						// cout << "Position Z: " << currentPosition[2] << endl << endl;

						// break;
					}
			}

			else if (pattern.substr(idx, 1).compare("+") == 0){													// ELSE IF the character matches +
				rotation = rand() % 65 + 1;																								// Chooses a random number to rotate the position left from 0 to 65.
				// rotation = 25;
				float numRotateZ =-  ((rotation * 3.14159) / 180);												// Variable - Converts degrees of the rotation to radians.
				rotateZ1[0] = cos(numRotateZ);																						// Updates the Z rotational matrix
				rotateZ1[1] = sin(numRotateZ);
				rotateZ1[4] = -sin(numRotateZ);
				rotateZ1[5] = cos(numRotateZ);
				multiply(rotateZ1, currentHeading, result1);															// Multiplies the rotateZ1 matrix by the currentHeading to get the new currentHeading.
				float magnitude = sqrt((result1[0]*result1[0]) + (result1[1]*result1[1]) + (result1[2]*result1[2]) + (result1[3]*result1[3]));	//Finds magnitude for normalization of heading.
					currentHeading[0] = result1[0] / magnitude;																// Normalizes the currentHeading vector.
					currentHeading[1] = result1[1] / magnitude;
					currentHeading[2] = result1[2] / magnitude;
					currentHeading[3] = result1[3] / magnitude;
				}

			else if (pattern.substr(idx, 1).compare("-") == 0){													// ELSE IF the character matches -
				rotation = rand() % 65 + 1;																								// Chooses a random number to rotate the position right from 0 to 65.
				// rotation = 25;
				float numRotateZ =+ ((rotation * 3.14159) / 180);													// Variable - Converts degrees of the rotation to radians.
				rotateZ1[0] = cos(numRotateZ);																						// Updates the Z rotational matrix
				rotateZ1[1] = sin(numRotateZ);
				rotateZ1[4] = -sin(numRotateZ);
				rotateZ1[5] = cos(numRotateZ);
				multiply(rotateZ1, currentHeading, result1);															// Multiplies the rotateZ1 matrix by the currentHeading to get the new currentHeading.
				float magnitude = sqrt((result1[0]*result1[0]) + (result1[1]*result1[1]) + (result1[2]*result1[2]) + (result1[3]*result1[3]));	//Finds magnitude for normalization of heading.
					currentHeading[0] = result1[0] / magnitude;																// Normalizes the currentHeading vector.
					currentHeading[1] = result1[1] / magnitude;
					currentHeading[2] = result1[2] / magnitude;
					currentHeading[3] = result1[3] / magnitude;
				}
			// }
		if (b1 == true) {

			// cout << "Points    " << pointsCount << endl;
			// cout << "Position X: " << currentPosition[0] << endl;
			// cout << "Position Y: " << currentPosition[1] << endl;
			// cout << "Position Z: " << currentPosition[2] << endl << endl;

			if (currentPosition[0] >= box || (currentPosition[1] >= box)){
				// cout << "HERE" << endl;

				lastPos[0] = currentPosition[0];
				lastPos[1] = currentPosition[1];
				lastPos[2] = currentPosition[2];
				lastPos[3] = currentPosition[3];

				cout << "Position X: " << currentPosition[0] << endl;
				cout << "Position Y: " << currentPosition[1] << endl;
				cout << "Position Z: " << currentPosition[2] << endl << endl;

				for (int times = 0; times < 3; times ++){
					currentPosition[0] -= (currentHeading[0]*.3);
					currentPosition[1] -= (currentHeading[1]*.3);
					currentPosition[2] -= (currentHeading[2]*.3);
				}


				// cout << "Before While" << endl;

				while (((currentPosition[0] < box) && (currentPosition[1] < box))){
					currentPosition[0] += currentHeading[0]*.01;
					currentPosition[1] += currentHeading[1]*.01;
					currentPosition[2] += currentHeading[2]*.01;

					// counter12++;
				}

				// cout << "After while" << endl;

				if ((currentPosition[0] < (box + 0.03)) && (currentPosition[1] < (box+0.03))){
					cout << "In if" << endl;
					branchPoints[pointsCount + 0] = currentPosition[0];												//Adds the currentPosition to the list of branching points
					branchPoints[pointsCount + 1] = currentPosition[1];
					branchPoints[pointsCount + 2] = currentPosition[2];
					pointsCount += 3;

					leafPoints[leafCount + 0] = currentPosition[0];														// Add the point to the leaf array (where leaves will be placed)
					leafPoints[leafCount + 1] = currentPosition[1];
					leafPoints[leafCount + 2] = currentPosition[2];
					leafCount += 3;
				}


				cout << "A Position X: " << currentPosition[0] << endl;
				cout << "Position Y: " << currentPosition[1] << endl;
				cout << "Position Z: " << currentPosition[2] << endl << endl;

				currentPosition[0] = lastPos[0];
				currentPosition[1] = lastPos[1];
				currentPosition[2] = lastPos[2];
				currentPosition[3] = lastPos[3];
			}

			}

			// cout << "HERE" << endl;
		}
	// }


	// cout << pointsCount << endl;
	// cout << totalCount << endl;
	GLfloat newPoints[pointsCount];

	for (int y = 0; y < pointsCount; y+=3){
		newPoints[y + 0] = branchPoints[y + 0];
		newPoints[y + 1] = branchPoints[y + 1];
		newPoints[y + 2] = branchPoints[y + 2];
		// cout << newPoints[y] << endl;
		// cout << newPoints[y + 1] << endl;
		// cout << newPoints[y + 2] << endl << endl;
	}


	// for (int y = 0; y < pointsCount; y += 3){
	// 	if( branchPoints[y] > 0.5){
	// 		cout << "Position X: " << branchPoints[y] << endl;
	// 		cout << "Position X: " << branchPoints[y + 1] << endl;
	// 		cout << "Position X: " << branchPoints[y + 2] << endl << endl;
	// 	}
	// 	else if (branchPoints[y + 1] > 0.5){
	// 		cout << "Position X: " << branchPoints[y] << endl;
	// 		cout << "Position X: " << branchPoints[y + 1] << endl;
	// 		cout << "Position X: " << branchPoints[y + 2] << endl << endl;
	// 	}
	// 	else if (branchPoints[y + 2] > 0.5){
	// 		cout << "Position X: " << branchPoints[y] << endl;
	// 		cout << "Position X: " << branchPoints[y + 1] << endl;
	// 		cout << "Position X: " << branchPoints[y + 2] << endl << endl;
	// 	}
	// }

	string modelName = "Leaf.obj";																								// Variable - Name of the OBJ to load.

	int numVert = countLabel(modelName, "v");																			// Variable - Counts the number of vertices.
	GLfloat* verts = new GLfloat[3*numVert];																			// Variable - creates array for the vertices
	loadVertices(modelName, verts);																								// Loads vertices into array

	int numFaces = countLabel(modelName, "f");																		// Variable - Counts the number of faces.
	GLint* faces = new GLint[3*numFaces];																					// Variable - creates array for faces
	loadFaces(modelName, faces);																									// Loads faces into array

	GLfloat* faceNormals = new GLfloat[3*numFaces];																// Variable - Creates array for face normals
	computeFaceNormals(faceNormals, verts, faces, numFaces);											// Computes the face normals in OBJ and adds them to the array

	GLfloat* vertNormals = new GLfloat[3*numVert];																// Variable - Creates array for vertex normals
	computeVertNormals(vertNormals, verts, numVert, faces, numFaces, faceNormals);// Computes the vertex normals in OBJ and adds them to the array

	int leavesWanted = leafCount / 3;																							// Variable - Number of leaves wanted for the tree
	GLfloat* points = new GLfloat[leavesWanted*9*numFaces];												// Variable - Creates array for the points of the OBJ
	GLfloat* normals = new GLfloat[leavesWanted*9*numFaces];											// Variable - Creates array for the normals of the OBJ

	for (int l = 0; l < leavesWanted; l++){																				// For loop - Makes multiple leaves based on leavesWanted
		for (int i = 0; i < numFaces; i++){																					// For loop - for all of the faces,
			int idx1 = faces[3*i + 0];																								// Variable - First index
			int idx2 = faces[3*i + 1];																								// Variable - Second index
			int idx3 = faces[3*i + 2];																								// Variable - Third index
			points[i*9 + 0 + l*29952] = verts[3*idx1+0];															// Continually adds the same numFaces values to fill points
			points[i*9 + 1 + l*29952] = verts[3*idx1+1];
			points[i*9 + 2 + l*29952] = verts[3*idx1+2];
			points[i*9 + 3 + l*29952] = verts[3*idx2+0];
			points[i*9 + 4 + l*29952] = verts[3*idx2+1];
			points[i*9 + 5 + l*29952] = verts[3*idx2+2];
			points[i*9 + 6 + l*29952] = verts[3*idx3+0];
			points[i*9 + 7 + l*29952] = verts[3*idx3+1];
			points[i*9 + 8 + l*29952] = verts[3*idx3+2];
			normals[i*9 + 0 + l*29952] = vertNormals[3*idx1+0];												// Continually adds the same numFaces values to fill normals
			normals[i*9 + 1 + l*29952] = vertNormals[3*idx1+1];
			normals[i*9 + 2 + l*29952] = vertNormals[3*idx1+2];
			normals[i*9 + 3 + l*29952] = vertNormals[3*idx2+0];
			normals[i*9 + 4 + l*29952] = vertNormals[3*idx2+1];
			normals[i*9 + 5 + l*29952] = vertNormals[3*idx2+2];
			normals[i*9 + 6 + l*29952] = vertNormals[3*idx3+0];
			normals[i*9 + 7 + l*29952] = vertNormals[3*idx3+1];
			normals[i*9 + 8 + l*29952] = vertNormals[3*idx3+2];
		}
	}
	int numPoints = 3*numFaces;																										// Variable - the number of points

	/*cout << "Leaf Count: " << leafCount << endl;
	cout << "numFaces: " << numFaces << endl << "numFaces-1: " << numFaces - 1 << endl;
	cout << "Total leaf array num: " << leavesWanted*9*numFaces << endl;
	cout << "Total XYZ for leaf points: " << leafCount << endl << "Total leaves: " << leafCount / 3 << endl << endl;*/

	cout << "\nCreating " << leafCount/3 << " leaves" << endl;										// Shows the user that it is creating leaves.
	for (int beginLeaf = 0; beginLeaf < leavesWanted; beginLeaf++) {							// For loop - Begins making multiple leaves.
																	 																							// Variable - Sets beginLeaf to 0 and counts up to the number of leaves needed.
		int endLeaf = beginLeaf + 1;																								// Variable - Sets endLeaf to be one greater than beginLeaf.
		//cout << "leaf " << endLeaf << endl;
		//cout << beginLeaf << " " << endLeaf << endl;															// Creates a "leaf" from beginLeaf to endLeaf.

		for (int i = beginLeaf*9*numFaces; i < endLeaf*9*numFaces - 1; i += 3){			// For loop - Starts loop to multiply each point by my matrices.
			int j;																																		// Variable - ??
			int k;																																		// Variable - ??

			float leafRotationZP;																											// Variable - Rotating leaf positively
			float leafRotationZN;																											// Variable - Rotating leaf negatively
			float leafRotationRadZP;																									// Variable - rotating leaf positively radians
			float leafRotationRadZN;																									// Variable - rotating leaf negatively radians

			GLfloat* new4x4 = new float[16];																					// Variable - Creates a new array to hold a 4x4 matrix.
			for (j = 0, k = 0; j < 16; j++){																					// For loop - Sets all values of 4x4 to 0.
				new4x4[j] = k;
			}

			GLfloat* newest4x4 = new float[16];																				// Variable - Creates a new array to hold a 4x4 matrix.
			for (j = 0, k = 0; j < 16; j++){																					// For loop - Sets all values of 4x4 to 0.
				newest4x4[j] = k;
			}

			GLfloat* newest4x4H = new float[16];																			// Variable - Creates a new array to hold a 4x4 matrix.
			for (j = 0, k = 0; j < 16; j++){																					// For loop - Sets all values of 4x4 to 0.
				newest4x4H[j] = k;
			}

			GLfloat* new4x1 = new float[4];																						// Variable - Creates a new array to hold a 4x1 vector.
			for (j = 0, k = 0; j < 4; j++){																						// For loop - Sets all values of 4x1 to 0.
				new4x1[j] = k;
			}

			dx = leafPoints[beginLeaf*3 + 0];																					// Gets the x value of the end of the current branch.
			dy = leafPoints[beginLeaf*3 + 1];																					// Gets the y value of the end of the current branch.
			dz = leafPoints[beginLeaf*3 + 2];																					// Gets the z value of the end of the current branch.

			float leafRot;																														// Variable - positive random number for rotation
			leafRot = RandomFloat(.1, .5);																						// Generates positive random number for leaf rotation
			// cout << leafRot << endl;

			float negLeafRot;																													// Variable - negative random number for rotation
			negLeafRot = RandomFloat(-.1, -.5);																				// Generates negative random number for leaf rotation
			// cout << negLeafRot << endl;

			leafRotationZP = leafRot*(1 - ((.5 - 0)/(1-0)) + 180*((.5 - 0)/(1-0)));		// ??
			leafRotationRadZP = -(leafRotationZP * 3.14159) / 180;										// Converts number to radians

			leafRotationZN = negLeafRot*(1 - ((.5 - 0)/(1-0)) + 180*((.5 - 0)/(1-0)));// ??
			leafRotationRadZN = -(leafRotationZN * 3.14159) / 180;										// Converts number to radians

			if (dx > 0) {																															// IF dx is greater than 0
				rotateZ3[0] = cos(leafRotationRadZP);																		// Set the third Z rotation matrix to be the positive rotation number
				rotateZ3[1] = sin(leafRotationRadZP);
				rotateZ3[4] = -sin(leafRotationRadZP);
				rotateZ3[5] = cos(leafRotationRadZP);
				// cout << "IF" << endl;
			} else if (dx == 0){																											// ELSE IF dx is equal to 0
				rotateZ3[0] = cos(0);																										// Set the third Z rotation matrix to be 0
				rotateZ3[1] = sin(0);
				rotateZ3[4] = -sin(0);
				rotateZ3[5] = cos(0);
				// cout << "else IF" << endl;
			} else {																																	// ELSE
				rotateZ3[0] = cos(leafRotationRadZN);																		// Set the third Z rotation matrix to be the negative rotation number
				rotateZ3[1] = sin(leafRotationRadZN);
				rotateZ3[4] = -sin(leafRotationRadZN);
				rotateZ3[5] = cos(leafRotationRadZN);
				// cout << "ELSE" <<endl;
			}

			translateMat[12] = dx;																										// Sets the dx value of translateMat to be the x-val of the current branch.
			translateMat[13] = dy;																										// Sets the dy value of translateMat to be the y-val of the current branch.
			translateMat[14] = dz;																										// Sets the dz value of translateMat to be the z-val of the current branch.

			scaleXNum = .0389;																												// Changes the scale of the leaf for x
			scaleYNum = .0389;																												// Changes the scale of the leaf for y
			scaleZNum = .0389;																												// Changes the scale of the leaf for z

			scale1[0] = scaleXNum;																										// Updates the scaling matrix
			scale1[5] = scaleYNum;																										// Updates the scaling matrix
			scale1[10] = scaleZNum;																										// Updates the scaling matrix

			numRotateZ2 = (25 * 3.14159) / -180;																			// Converts 25 degrees to radians

			rotateZ2[0] = cos(numRotateZ2);																						// Updates the second Z rotational matrix with 25 degrees to radians
			rotateZ2[1] = sin(numRotateZ2);
			rotateZ2[5] = -sin(numRotateZ2);
			rotateZ2[6] = cos(numRotateZ2);

			new4x4 = multiplyAgain(rotateX1, scale1, new4x4);													// Multiplies two 4x4 matrices together and makes a new 4x4 matrix.
			newest4x4H = multiplyAgain(rotateZ3, new4x4, newest4x4H);									// Multiplies two 4x4 matrices together and makes a new 4x4 matrix
			newest4x4 = multiplyAgain(translateMat, newest4x4H, newest4x4);						// Multiplies two 4x4 matrices together and makes a new 4x4 matrix.

			float currentLeafPoint[] = {points[i], points[i + 1], points[i + 2], 1};	// Variable - Gets the x, y, z values from points (leaf) to multiply by.

			multiply(newest4x4, currentLeafPoint, new4x1);														// Multiplies a 4x4 and a 4x1 matrix together and makes a new 4x1 matrix.

			points[i + 0] = new4x1[0];																								// Sets current points x-val to be the x-val of the new	4x1 matrix.
			points[i + 1] = new4x1[1];																								// Sets current points y-val to be the y-val of the new	4x1 matrix.
			points[i + 2] = new4x1[2];																								// Sets current points z-val to be the z-val of the new	4x1 matrix.
		}
		//cout << "Rotation Amount " << leafRotationRad << endl;
		//cout << dx << " " << dy << " " << dz << endl;
		//cout << beginLeaf*9*numFaces << " " << endLeaf*9*numFaces - 1 << endl << endl;
	}
	cout << "Done creating " << leafCount / 3 << " leaves\n" << endl;							// Prints out the number of leaves that are being created.

	/* these are the strings of code for the shaders
	the vertex shader positions each vertex point */
	const char *vertex_shader = "#version 410\n"																	// Vertex Shader for tree
		"out vec3 colour;"
		// "attribute vec3 vp;"
		"layout (location = 0) in vec3 vp;"
		"layout (location = 1) in vec3 vertex_normal;"
		"uniform mat4 model, view, proj;"
		"out vec3 position_eye, normal_eye;"
		"void main () {"
		"	position_eye = vp;"
		"	normal_eye = vertex_normal, 1.0;"
		"  gl_Position = proj * view * model * vec4(vp, 1.0);"											// ADDING IN *ORTHO	BREAKS THE LEAVES FROM THE BRANCHES	// Tree position.
		"  colour = vec3 (255, 0, 0);"
		"}";


	const char *geometry_shader = "#version 410\n"																// Geometry Shader for tree - creates the thicker branches
		//	MAKE SURE 'in' are arrays
		"layout (lines) in;" // lines, line_Strip is output, not input.
		// convert to points, line_strip, or triangle_strip
		"layout (triangle_strip, max_vertices = 6) out;"
		// NB: in and out pass-through vertex->fragment variables must go here if used
		"in vec3 colour[];"
		"out vec3 f_colour;"
		"out vec3 position_eye;"
		"out vec3 normal_eye;"

		"void thickLine(vec4 position, vec4 position2, float top, float bottom){"
			"gl_Position = position + vec4(-bottom, -bottom, 0.0, 0.0);" 							//bottom left
			"EmitVertex();"
			"gl_Position = position + vec4(bottom, -bottom, 0.0, 0.0);" 							//bottom right
			"EmitVertex();"
			"gl_Position = position2 + vec4(-top, top, 0.0, 0.0);" 										//top left
			"EmitVertex();"
			"gl_Position = position2 + vec4(top, top, 0.0, 0.0);" 										//top right
			"EmitVertex();"
			"EndPrimitive();"
		"}"

		"void main () {"
			"float bottom = 0.009;"
			"float top = 0.007;"
			"for(int i = 0; i < gl_in.length(); i++) {"
						"thickLine(gl_in[0].gl_Position, gl_in[1].gl_Position, top, bottom);"
			"}"
		"}";

	/* the fragment shader colours each fragment (pixel-sized area of the
	triangle) */
	const char *fragment_shader = "#version 410\n"																// Fragment Shader for tree - adds color and shading to the tree branches
		"in vec3 position_eye, normal_eye;"																					// Light comes from using the Phong model
		// fixed point light properties
		"vec3 light_position_world  = vec3 (2.0, 1.0, 0.0);"
		"vec3 Ls = vec3 (1.0, 0.0, 0.0);" // white specular colour
		"vec3 Ld = vec3 (0.7, 0.7, 0.7);" // dull white diffuse light colour
		"vec3 La = vec3 (0.2, 0.2, 0.2);" // grey ambient colour

		// surface reflectance
		"vec3 Ks = vec3 (1.0, 0.0, 0.0);" // fully reflect specular light
		"vec3 Kd = vec3 (0.0, 1.0, 0.0);" // green diffuse surface reflectance
		"vec3 Ka = vec3 (1.0, 1.0, 1.0);" // fully reflect ambient light
		"float specular_exponent = 100.0;" // specular 'power'

		"in vec3 f_colour;"
		"out vec4 frag_colour;"
		"void main () {"
		// ambient intensity
		"vec3 Ia = La * Ka;"

		// diffuse intensity
		// raise light position to eye space
		"vec3 light_position_eye = vec3 (vec4 (light_position_world, 1.0));"
		"vec3 distance_to_light_eye = light_position_eye - position_eye;"
		"vec3 direction_to_light_eye = normalize (distance_to_light_eye);"
		"float dot_prod = dot (direction_to_light_eye, normal_eye);"
		"dot_prod = max (dot_prod, 0.0);"
		"vec3 Id = Ld * Kd * dot_prod;" // final diffuse intensity

		// specular intensity
		"vec3 surface_to_viewer_eye = normalize (-position_eye);"

		// blinn
		"vec3 half_way_eye = normalize (surface_to_viewer_eye + direction_to_light_eye);"
		"float dot_prod_specular = max (dot (half_way_eye, normal_eye), 0.0);"
		"float specular_factor = pow (dot_prod_specular, specular_exponent);"

		"vec3 Is = Ls * Ks * specular_factor;" // final specular intensity

		"  frag_colour = vec4 (Ia + Id + Is, 1.0);"
		"  frag_colour = vec4 (0.545, 0.27, 0.074, 1.0);"    												//Makes tree brown.
		"}";

	/* GL shader objects for vertex and fragment shader [components] */
	GLuint vert_shader, frag_shader, geoShader;																		// Variables - creates three shader variables
	/* GL shader program object [combined, to link] */
	GLuint shader_programme;																											// Variable - creates shader program variable

	//---------------------------------------------------------------------------------------------------------------------------  SHADERS FOR THE LEAF
	/* these are the strings of code for the shaders
	the vertex shader positions each vertex point */
	const char *vertex_shader2 = "#version 410\n"																	// Vertex Shader for leaf.
	//	"attribute vec3 vp;"
		"layout (location = 0) in vec3 vp;"
		"layout (location = 1) in vec3 vertex_normal;"
	  "uniform mat4 model, view, proj;" 																					// Gets matrices inside vertex shader.
		"out vec3 position_eye, normal_eye;"
		"void main () {"
		"	position_eye = vp;"
		"	normal_eye = vertex_normal, 1.0;"																					// Multiplies vec4 by matrices.
		"  gl_Position = proj * view * model * vec4(vp, 1.0);"											// ADDING IN *ORTHO	BREAKS THE LEAVES FROM THE BRANCHES	// Tree position.
		"}";

	/* the fragment shader colours each fragment (pixel-sized area of the
	triangle) */
	const char *fragment_shader2 = "#version 410\n"																// Fragment Shader - adds in lighting and color of leaves using the Phong model
		"in vec3 position_eye, normal_eye;"

		// fixed point light properties
		"vec3 light_position_world  = vec3 (2.0, 1.0, 0.0);"
		"vec3 Ls = vec3 (1.0, 0.0, 0.0);" // white specular colour
		"vec3 Ld = vec3 (0.7, 0.7, 0.7);" // dull white diffuse light colour
		"vec3 La = vec3 (0.2, 0.2, 0.2);" // grey ambient colour

		// surface reflectance
		"vec3 Ks = vec3 (1.0, 0.0, 0.0);" // fully reflect specular light
		"vec3 Kd = vec3 (0.0, 1.0, 0.0);" // green diffuse surface reflectance
		"vec3 Ka = vec3 (1.0, 1.0, 1.0);" // fully reflect ambient light
		"float specular_exponent = 100.0;" // specular 'power'

		"out vec4 frag_colour;"
		"void main () {"
		// ambient intensity
		"vec3 Ia = La * Ka;"

		// diffuse intensity
		// raise light position to eye space
		"vec3 light_position_eye = vec3 (vec4 (light_position_world, 1.0));"
		"vec3 distance_to_light_eye = light_position_eye - position_eye;"
		"vec3 direction_to_light_eye = normalize (distance_to_light_eye);"
		"float dot_prod = dot (direction_to_light_eye, normal_eye);"
		"dot_prod = max (dot_prod, 0.0);"
		"vec3 Id = Ld * Kd * dot_prod;" // final diffuse intensity

		// specular intensity
		"vec3 surface_to_viewer_eye = normalize (-position_eye);"

		// blinn
		"vec3 half_way_eye = normalize (surface_to_viewer_eye + direction_to_light_eye);"
		"float dot_prod_specular = max (dot (half_way_eye, normal_eye), 0.0);"
		"float specular_factor = pow (dot_prod_specular, specular_exponent);"

		"vec3 Is = Ls * Ks * specular_factor;" // final specular intensity

		"  frag_colour = vec4 (Ia + Id + Is, 1.0);"
		"}";
	/* GL shader objects for vertex and fragment shader [components] */
	GLuint vert_shader2, frag_shader2;																						// Variables - creates shader variables
	/* GL shader program object [combined, to link] */
	GLuint shader_programme2;																											// Variable - creates second shader program variable

	//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	/* start GL context and O/S window using the GLFW helper library */
	if ( !glfwInit() ) {
		fprintf( stderr, "ERROR: could not start GLFW3\n" );
		return 1;
	}

	g_window = glfwCreateWindow(640, 480, "Fractal Window", NULL, NULL);
	if ( !g_window ) {
		fprintf( stderr, "ERROR: could not open window with GLFW3\n" );
		glfwTerminate();
		return 1;
	}

	glfwMakeContextCurrent( g_window );
	/* start GLEW extension handler */
	glewExperimental = GL_TRUE;
	glewInit();

	/* get version info */
	renderer = glGetString( GL_RENDERER ); /* get renderer string */
	version = glGetString( GL_VERSION );	 /* version as a string */
	printf( "Renderer: %s\n", renderer );
	printf( "OpenGL version supported %s\n", version );

	/* tell GL to only draw onto a pixel if the shape is closer to the viewer
	than anything already drawn at that pixel */
	glEnable( GL_DEPTH_TEST ); /* enable depth-testing */
	/* with LESS depth-testing interprets a smaller depth value as meaning "closer" */
	glDepthFunc( GL_LESS );
	/* a vertex buffer object (VBO) is created here. this stores an array of
	data on the graphics adapter's memory. in our case - the vertex points */
	glGenBuffers( 1, &vbo );
	glBindBuffer( GL_ARRAY_BUFFER, vbo );
	glBufferData( GL_ARRAY_BUFFER, (pointsCount) * sizeof( GLfloat ), newPoints, GL_STATIC_DRAW ); // count*3 to not lose points

	/* the vertex array object (VAO) is a little descriptor that defines which
	data from vertex buffer objects should be used as input variables to vertex
	shaders. in our case - use our only VBO, and say 'every three floats is a
	variable' */
	glGenVertexArrays( 1, &vao );
	glBindVertexArray(vao);
	// "attribute #0 should be enabled when this vao is bound"
	glEnableVertexAttribArray(0);
	// this VBO is already bound, but it's a good habit to explicitly specify which
	// VBO's data the following
	// vertex attribute pointer refers to
	glBindBuffer( GL_ARRAY_BUFFER, vbo );
	// "attribute #0 is created from every 3 variables in the above buffer, of type
	// float (i.e. make me vec3s)"
	glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, NULL );

	//--------------------------------------------------------------------------------------------------------------------------------------------  LEAF SHADERS
	GLuint points_vbo;
	glGenBuffers (1, &points_vbo);
	glBindBuffer (GL_ARRAY_BUFFER, points_vbo);
	glBufferData (GL_ARRAY_BUFFER, leavesWanted * 3 * numPoints * sizeof (GLfloat), points, GL_STATIC_DRAW);

	GLuint normals_vbo;
	glGenBuffers (1, &normals_vbo);
	glBindBuffer (GL_ARRAY_BUFFER, normals_vbo);
	glBufferData (GL_ARRAY_BUFFER, leavesWanted * 3 * numPoints * sizeof (GLfloat), normals, GL_STATIC_DRAW);

	GLuint vao2;
	glGenVertexArrays (1, &vao2);
	glBindVertexArray (vao2);
	glBindBuffer (GL_ARRAY_BUFFER, points_vbo);
	glVertexAttribPointer (0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
	glBindBuffer (GL_ARRAY_BUFFER, normals_vbo);
	glVertexAttribPointer (1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
	glEnableVertexAttribArray (0);
	glEnableVertexAttribArray (1);

	//--------------------------------------------------------------------------------------------------------------------------------------------

	/* here we copy the shader strings into GL shaders, and compile them. we
	then create an executable shader 'program' and attach both of the compiled
	shaders. we link this, which matches the outputs of the vertex shader to
	the inputs of the fragment shader, etc. and it is then ready to use */
	vert_shader = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vert_shader, 1, &vertex_shader, NULL);
	glCompileShader(vert_shader);

	int params = -1;
		glGetShaderiv( vert_shader, GL_COMPILE_STATUS, &params );
		if ( GL_TRUE != params ) {
			fprintf( stderr, "ERROR: vert GL shader index %i did not compile\n", vert_shader );
			return 1; // or exit or something
		}

	geoShader = glCreateShader(GL_GEOMETRY_SHADER);
	glShaderSource(geoShader, 1, &geometry_shader, NULL);
	glCompileShader(geoShader);

	// check for compile errors
		glGetShaderiv( geoShader, GL_COMPILE_STATUS, &params );
		if ( GL_TRUE != params ) {
			fprintf( stderr, "ERROR: geom GL shader index %i did not compile\n", geoShader );
			return 1; // or exit or something
		}

	frag_shader = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(frag_shader, 1, &fragment_shader, NULL);
	glCompileShader(frag_shader);

	// check for compile errors
		glGetShaderiv( frag_shader, GL_COMPILE_STATUS, &params );
		if ( GL_TRUE != params ) {
			fprintf( stderr, "ERROR: frag GL shader index %i did not compile\n", frag_shader );
			return 1; // or exit or something
		}

	shader_programme = glCreateProgram();
	glAttachShader(shader_programme, frag_shader);
	glAttachShader(shader_programme, geoShader);
	glAttachShader(shader_programme, vert_shader);
	glLinkProgram(shader_programme);
	// glPointSize(5.0);

	glGetProgramiv( shader_programme, GL_LINK_STATUS, &params );
	if ( GL_TRUE != params ) {
		fprintf( stderr, "ERROR: could not link shader_programme GL index %i\n",
						 shader_programme );
		return false;
	}
	( is_programme_valid( shader_programme ) );

	//-----------------------------------------------------------------------------------------------------------------------------------------  Leaf Shaders
	vert_shader2 = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vert_shader2, 1, &vertex_shader2, NULL);
	glCompileShader(vert_shader2);

	// int params = -1;
		glGetShaderiv( vert_shader2, GL_COMPILE_STATUS, &params );
		if ( GL_TRUE != params ) {
			fprintf( stderr, "ERROR: vert GL shader index %i did not compile\n", vert_shader2 );
			return 1; // or exit or something
		}

	frag_shader2 = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(frag_shader2, 1, &fragment_shader2, NULL);
	glCompileShader(frag_shader2);

	// check for compile errors
		glGetShaderiv( frag_shader2, GL_COMPILE_STATUS, &params );
		if ( GL_TRUE != params ) {
			fprintf( stderr, "ERROR: frag GL shader index %i did not compile\n", frag_shader2 );
			return 1; // or exit or something
		}


	shader_programme2 = glCreateProgram();
	glAttachShader(shader_programme2, frag_shader2);
	glAttachShader(shader_programme2, vert_shader2);
	glLinkProgram(shader_programme2);

	glGetProgramiv( shader_programme2, GL_LINK_STATUS, &params );
	if ( GL_TRUE != params ) {
		fprintf( stderr, "ERROR: could not link shader_programme GL index %i\n",
						 shader_programme2 );
		return false;
	}
	( is_programme_valid( shader_programme2 ) );


	//--------------------------------------------------------------------------------------------------------------------------------------------------------
	/* this loop clears the drawing surface, then draws the geometry described
			by the VAO onto the drawing surface. we 'poll events' to see if the window
			was closed, etc. finally, we 'swap the buffers' which displays our drawing
			surface onto the view area. we use a double-buffering system which means
			that we have a 'currently displayed' surface, and 'currently being drawn'
			surface. hence the 'swap' idea. in a single-buffering system we would see
			stuff being drawn one-after-the-other */
  multiplyNew(view, proMat, viewResult);
	while ( !glfwWindowShouldClose( g_window ) ) {
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		multiplyNew(rotateY, rotateX, resultRotation);
		multiplyNew(rotateZ, resultRotation, result);
		multiplyNew(result, trans, result);
		multiplyNew(result, translate, result);

		//View matrix info
		int trans_mat_location2 = glGetUniformLocation (shader_programme2, "model");
		glUseProgram( shader_programme2 );
		glUniformMatrix4fv (trans_mat_location2, 1, GL_FALSE, result);
		int view_mat_location2 = glGetUniformLocation (shader_programme2, "view");
		glUseProgram( shader_programme2 );
		glUniformMatrix4fv (view_mat_location2, 1, GL_FALSE, lookAt);
		int proj_mat_location2 = glGetUniformLocation (shader_programme2, "proj");
		glUseProgram( shader_programme2 );
		glUniformMatrix4fv (proj_mat_location2, 1, GL_FALSE, viewResult);
		//glUniformMatrix4fv (proj_mat_location, 1, GL_FALSE, proMat);


		glUseProgram(shader_programme2);
		glBindVertexArray(vao2);
		// glDrawArrays(GL_TRIANGLES, 0, leavesWanted * numPoints);

		//View matrix info
		int trans_mat_location = glGetUniformLocation (shader_programme, "model");
		glUseProgram( shader_programme );
		glUniformMatrix4fv (trans_mat_location, 1, GL_FALSE, result);
		int view_mat_location = glGetUniformLocation (shader_programme, "view");
		glUseProgram( shader_programme );
		glUniformMatrix4fv (view_mat_location, 1, GL_FALSE, lookAt);
		int proj_mat_location = glGetUniformLocation (shader_programme, "proj");
		glUseProgram( shader_programme );
		glUniformMatrix4fv (proj_mat_location, 1, GL_FALSE, viewResult);
		//glUniformMatrix4fv (proj_mat_location, 1, GL_FALSE, proMat);

		/* wipe the drawing surface clear */
		glUseProgram(shader_programme);
		glBindVertexArray(vao);

		/* draw points 0-3 from the currently bound VAO with current in-use shader */
		glDrawArrays(GL_LINE_STRIP, 0, totalCount); //GL_LINE_STRIP

		/* update other events like input handling */
		glfwPollEvents();

		glfwSetKeyCallback(g_window, key_callback);

		translate[12]=x;
		translate[13]=y;
		translate[14]=z;
		scale[0]=sx;
		scale[5]=sy;
		scale[10]=sz;
		rotateX[5]=cos(rx);
		rotateX[6]=sin(rx);
		rotateX[9]=-sin(rx);
		rotateX[10]=cos(rx);
		rotateY[0]=cos(ry);
		rotateY[2]=-sin(ry);
		rotateY[8]=sin(ry);
		rotateY[10]=cos(ry);
		rotateZ[0]=cos(rz);
		rotateZ[1]=sin(rz);
		rotateZ[4]=-sin(rz);
		rotateZ[5]=cos(rz);
		float range = tan(fov*0.5)*near;
		proMat[0] = (2*near)/((range*aspect)+(range*aspect));
		proMat[5] = near/range;
		proMat[10] = -(far+near)/(far-near);
		proMat[14] = -(2*far*near)/(far-near);
		viewMat = ortho;
		multiplyNew(view, viewMat, viewResult);

		/* put the stuff we've been drawing onto the display */
		glfwSwapBuffers( g_window );
	}

	/* close GL context and any other GLFW resources */
	glfwTerminate();
	return 0;
}
