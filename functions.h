#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_
#include <string>
#include <string.h>
void multiply(GLfloat, GLfloat, GLfloat);

GLfloat* multiplyAgain(GLfloat, GLfloat, GLfloat);

string generatePattern();

int countF(string);

int countbracket(string);

int countLabel(string, char);

void loadVertices(string, GLfloat);

void computeFaceNormals(GLfloat, GLfloat, GLfloat, int);

void computeVertNormals(GLfloat, GLfloat, int, GLint, int, GLfloat);

void loadFaces(string, GLint);

#endif
