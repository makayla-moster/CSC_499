//g++ -w -o Makayla.exe gl_utils.cpp maths_funcs.cpp FractalPieces.cpp libglfw3dll.a libglew32.dll.a -I include -lglfw3 -lgdi32 -lopengl32

//g++ -w -o Makayla.exe FractalPieces.cpp libglfw3dll.a libglew32.dll.a -I include -lglfw3 -lgdi32 -lopengl32

#include "gl_utils.h"
#include "maths_funcs.h"
#include <GL/glew.h>		/* include GLEW and new version of GL on Windows */
#include <GLFW/glfw3.h>     /* GLFW helper library */
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <stack>
#include <vector>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <time.h>
#include <assert.h>
#include <string>
#include <string.h>
#include <stdarg.h> 
#define GL_LOG_FILE "gl.log"  
using namespace std;

/*void multiply(GLfloat matrix1[], GLfloat matrix2[], GLfloat result[]){   // this is to multiply a 4x4 and a 4x1 matrix together
	result[0] = (matrix1[0]*matrix2[0]) + (matrix1[4]*matrix2[1]) + (matrix1[8]*matrix2[2]) + (matrix1[12]*matrix2[3]);
	result[1] = (matrix1[1]*matrix2[0]) + (matrix1[5]*matrix2[1]) + (matrix1[9]*matrix2[2]) + (matrix1[13]*matrix2[3]);
	result[2] = (matrix1[2]*matrix2[0]) + (matrix1[6]*matrix2[1]) + (matrix1[10]*matrix2[2]) + (matrix1[14]*matrix2[3]);
	result[3] = (matrix1[3]*matrix2[0]) + (matrix1[7]*matrix2[1]) + (matrix1[11]*matrix2[2]) + (matrix1[15]*matrix2[3]);
}*/

/*GLfloat* multiplyAgain(GLfloat matrix1[], GLfloat matrix2[], GLfloat result[]){
	for (int i = 0; i < 16; i++){
		for (int j = 0; j < 4; j++){
			result[i] = result[i] + (matrix2[(i / 4) * 4 + j] * matrix1[(i % 4)+(4*j)]);
			//cout << matrix2[(i / 4) * 4 + j]  * matrix1[(i % 4)+(4*j)] << endl;
			//cout << i << " " << j << " " << result[i] << endl;
		}
	}
	return result;
}*/

string generatePattern(){												//Generates a pattern to create a tree.
    int numIts = 1; // Number of iterations
    string pattern = "F"; //"[X]";    // Using F for the pattern 
    
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
	int count = 0;
	for (int idx = 0; idx < pattern.length(); idx++){					// Loops through string pattern to find number of 'F' within the pattern.
		if (pattern.substr(idx, 1).compare("F") == 0){
			count ++;
		}
	}
	return count;
}

int countbracket(string pattern) {
	int countBracket = 0;
	for (int idx = 0; idx < pattern.length(); idx++){					// Loops through string pattern to find number of ']' within the pattern.
		if (pattern.substr(idx, 1).compare("]") == 0){
			countBracket ++;
		}
	}	
	
	return countBracket;
}

int countLabel(string modelName, char label[]){
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

/*void loadVertices(string modelName, GLfloat verts[]){	
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
}*/


/*void computeFaceNormals(GLfloat faceNormals[], GLfloat verts[], GLint faces[], int numFaces){
	for (int i = 0; i < numFaces; i++){
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
	
}*/

/*void computeVertNormals(GLfloat normals[], GLfloat verts[], int numVerts, GLint faces[], int numFaces, GLfloat faceNormals[]){
	for (int i = 0; i < numVerts; i++){
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
}*/

/*void loadFaces(string modelName, GLint faces[]){    					//To read in Maya OBJ files.
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
}*/


int g_gl_width = 640;
int g_gl_height = 480;
GLFWwindow* g_window = NULL;

int main(){
	const GLubyte *renderer;
	const GLubyte *version;
	GLuint vao;
	GLuint vbo;
	
	
	GLfloat* resultAgain = new float[16];
	
	int rotation;													// Rotation of branches.
	int count;														// Counts number of 'F' in string.
	int countBracket;												// Counts number of ']' in string.
	
	string pattern = generatePattern();								// Generates string pattern to make tree from.
	cout << pattern << endl << endl;
	count = countF(pattern);										// Function to count the number of 'F' in the string.
	cout << count << endl;
	countBracket = countbracket(pattern);							// Function to count the number of ']' in the string.
	cout << countBracket << endl;
	int totalCount = count + countBracket;							//Total amount of points, including the backtracking points that are added for the lines.
	
	GLfloat branchPoints[totalCount*3]; 							//List of points to make the branches. Includes extra points for lines.
	GLfloat leafPoints[countBracket*3];								//List of points to place the leaves - only on the ends of branches though.
	GLfloat* result = new float[4];									//This is a new 4x1 matrix to acquire the new heading.	
	stack<float> PositionStack;										//Stack to put branch point positions on.
	stack<float> HeadingStack;										//Stack to put branch headings on.
	int pointsCount = 0;											//Counts number of points needed to make a matrix of points.
	int leafCount = 0;												//Counts number of points needed to make a matrix for leaves.
	
	float rz = (90 * 3.14159) / 180;								//For the first rotation, so the trunk is 90 degrees from the bottom of the screen. Converts degrees to radians.
	float rz2 = 0; //(30 * 3.14159) / 180;								//For rotation of leaf.
	float rx = -(90 * 3.14159) / 180;								//Rotates the leaf OBJ to be upright in radians.
	float ry = (0 * 3.14159) / 180;									//Rotates the leaf along the y-axis.
	GLfloat dx = 0;													//For leaf translation along the x-axis.
	GLfloat dy = 0;													//For leaf translation along the y-axis.
	GLfloat dz = 0;													//For leaf translation along the z-axis.
	GLfloat sx = 1; //.0095;										//Scales the leaf.
	GLfloat sy = 1; //.1095;										//Scales the leaf.
	GLfloat sz = 1; //.1095;										//Scales the leaf.
	GLfloat currentPosition[] = {0.0f, -0.25f, 0.0f, 1.0f};			//Beginning current position of the tree.
	GLfloat currentHeading[] = {0.0f, 0.5f, 0.0f, 0.0f};			//Beginning current heading of the tree.
	GLfloat rotateZ[] = 											//Rotation matrix for the z-axis.
		{cos(rz),sin(rz),0,0,
		-sin(rz),cos(rz),0,0,
		0,0,1,0,
		0,0,0,1};
	GLfloat rotateZ2[] = 																		//Rotation matrix for the z-axis.
		{cos(rz2),sin(rz2),0,0,
		-sin(rz2),cos(rz2),0,0,
		0,0,1,0,
		0,0,0,1};
	GLfloat scale[] =																			//Scale matrix.
		{sx,0,0,0,
		 0,sy,0,0,
		 0,0,sz,0,
		 0,0,0,1};
	GLfloat rotateX[] = 																		//Rotation matrix for the x-axis.
		{1,0,0,0,
		 0,cos(rx),sin(rx),0,
		 0,-sin(rx),cos(rx),0,
		 0,0,0,1}; 
	GLfloat rotateY[] =																			//Rotation matrix for the y-axis. 
		{cos(ry),0,-sin(ry),0,
		 0,1,0,0,
		 sin(ry),0,cos(ry),0,
		 0,0,0,1};
	GLfloat translateMat[] =																		//Translation matrix.
		{1,0,0,0,
		 0,1,0,0,
		 0,0,1,0,
		 dx,dy,dz,1};
	GLfloat identity[] = 
		{1,0,0,0,
		 0,1,0,0,
		 0,0,1,0,
		 0,0,0,1};
		 
	branchPoints[pointsCount + 0] = currentPosition[0];											//These lines add the first set of points to the list of points.
	branchPoints[pointsCount + 1] = currentPosition[1];
	branchPoints[pointsCount + 2] = currentPosition[2];
	pointsCount += 3;
	
	//cout << pointsCount << endl;
	
		/*for (int idx = 0; idx < pattern.length(); idx++){											//Parser to parse through the Tree String pattern.
			/*if (pattern.substr(idx,1).compare("[") == 0){
				PositionStack.push(currentPosition[3]);												//Pushes the currentPosition onto the PositionStack.
				PositionStack.push(currentPosition[2]);
				PositionStack.push(currentPosition[1]);
				PositionStack.push(currentPosition[0]);
				
				HeadingStack.push(currentHeading[3]);												//Pushes the currentHeading onto the HeadingStack.
				HeadingStack.push(currentHeading[2]);
				HeadingStack.push(currentHeading[1]);
				HeadingStack.push(currentHeading[0]);
			}*/
			
			/*else if (pattern.substr(idx, 1).compare("]") == 0){
				
				leafPoints[leafCount + 0] = currentPosition[0];										//Adds the point to the array which will be where I put my leaf OBJ.
				leafPoints[leafCount + 1] = currentPosition[1];
				leafPoints[leafCount + 2] = currentPosition[2];
				leafCount += 3;
				
				cout << endl;
				cout << "Position X: " << currentPosition[0] << endl;
				cout << "Position Y: " << currentPosition[1] << endl;
				cout << "Position Z: " << currentPosition[2] << endl << endl;
				
				currentPosition[0] = PositionStack.top();											//Sets the current position back to the top of the stack.
				PositionStack.pop();																//Pops the current position from the top of the stack.
				currentPosition[1] = PositionStack.top();
				PositionStack.pop();
				currentPosition[2] = PositionStack.top();
				PositionStack.pop();
				currentPosition[3] = PositionStack.top();
				PositionStack.pop();
				
				branchPoints[pointsCount + 0] = currentPosition[0];									//Adds the currentPosition to the list of points.
				branchPoints[pointsCount + 1] = currentPosition[1];
				branchPoints[pointsCount + 2] = currentPosition[2];
				pointsCount += 3;
			
				currentHeading[0] = HeadingStack.top();												//Sets the currentHeading to the top of the HeadingStack.
				HeadingStack.pop();																	//Pops the currentHeading from the top of the stack.
				currentHeading[1] = HeadingStack.top();
				HeadingStack.pop();
				currentHeading[2] = HeadingStack.top();
				HeadingStack.pop();
				currentHeading[3] = HeadingStack.top();
				HeadingStack.pop();
			}*/
			
			/*else if (pattern.substr(idx, 1).compare("F") == 0){			
				currentPosition[0] += currentHeading[0]*.2;											//Changes the height of the tree, I like .2.
				currentPosition[1] += currentHeading[1]*.2;
				currentPosition[2] += currentHeading[2]*.2;
				currentPosition[3] += currentHeading[3];
				
				branchPoints[pointsCount + 0] = currentPosition[0];									//Adds the currentPosition to the list of points.
				branchPoints[pointsCount + 1] = currentPosition[1];
				branchPoints[pointsCount + 2] = currentPosition[2];
				pointsCount += 3;
			}*/
			
			/*else if (pattern.substr(idx, 1).compare("+") == 0){
				rotation = rand() % 65 + 1;															//Chooses a random number to rotate by from 0 to 65.
				//rotation = 25;
				float rz =-  ((rotation * 3.14159) / 180);											//Converts degrees of the rotation to radians.
				rotateZ[0] = cos(rz);
				rotateZ[1] = sin(rz);
				rotateZ[4] = -sin(rz);
				rotateZ[5] = cos(rz);
				multiply(rotateZ, currentHeading, result);											//Multiplies the rotateZ matrix by the currentHeading to get the new currentHeading.
				float magnitude = sqrt((result[0]*result[0]) + (result[1]*result[1]) + (result[2]*result[2]) + (result[3]*result[3]));	//Finds magnitude for normalization of heading.
				currentHeading[0] = result[0] / magnitude;											//Normalizes the currentHeading vector.
				currentHeading[1] = result[1] / magnitude;
				currentHeading[2] = result[2] / magnitude;
				currentHeading[3] = result[3] / magnitude;
			}*/
			
			/*else if (pattern.substr(idx, 1).compare("-") == 0){
				rotation = rand() % 65 + 1;															//Chooses a random number to rotate by from 0 to 65.
				//rotation = 25;
				float rz =+ ((rotation * 3.14159) / 180);											//Converts degrees of the rotation to radians.
				rotateZ[0] = cos(rz);
				rotateZ[1] = sin(rz);
				rotateZ[4] = -sin(rz);
				rotateZ[5] = cos(rz);
				multiply(rotateZ, currentHeading, result);											//Multiplies the rotateZ matrix by the currentHeading to get the new currentHeading.
				float magnitude = sqrt((result[0]*result[0]) + (result[1]*result[1]) + (result[2]*result[2]) + (result[3]*result[3]));	//Finds magnitude for normalization of heading.
				currentHeading[0] = result[0] / magnitude;											//Normalizes the currentHeading vector.
				currentHeading[1] = result[1] / magnitude;
				currentHeading[2] = result[2] / magnitude;
				currentHeading[3] = result[3] / magnitude;
			}*/
		//}
	
	
	/* these are the strings of code for the shaders
	the vertex shader positions each vertex point */
	const char *vertex_shader = "#version 410\n"												// Vertex Shader for tree.
		"attribute vec3 vp;"
		"void main () {"
		"  gl_Position = vec4 (vp, 1.0);"														// Tree position.
		"}";
	/* the fragment shader colours each fragment (pixel-sized area of the
	triangle) */
	const char *fragment_shader = "#version 410\n"												// Fragment Shader for tree.
		"out vec4 frag_colour;"
		"void main () {"
		"  frag_colour = vec4 (0.545, 0.27, 0.074, 1.0);"    									//Makes tree brown.
		"}";
	/* GL shader objects for vertex and fragment shader [components] */
	GLuint vert_shader, frag_shader;
	/* GL shader program object [combined, to link] */
	GLuint shader_programme;
	
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
	glBufferData( GL_ARRAY_BUFFER, (totalCount*3) * sizeof( GLfloat ), branchPoints, GL_STATIC_DRAW ); // count*3 to not lose points 

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
	
	/* here we copy the shader strings into GL shaders, and compile them. we
	then create an executable shader 'program' and attach both of the compiled
	shaders. we link this, which matches the outputs of the vertex shader to
	the inputs of the fragment shader, etc. and it is then ready to use */
	vert_shader = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vert_shader, 1, &vertex_shader, NULL);
	glCompileShader(vert_shader);
	frag_shader = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(frag_shader, 1, &fragment_shader, NULL);
	glCompileShader(frag_shader);
	shader_programme = glCreateProgram();
	glAttachShader(shader_programme, frag_shader);
	glAttachShader(shader_programme, vert_shader);
	glLinkProgram(shader_programme);
	glPointSize(5.0);
	
	/* this loop clears the drawing surface, then draws the geometry described
			by the VAO onto the drawing surface. we 'poll events' to see if the window
			was closed, etc. finally, we 'swap the buffers' which displays our drawing
			surface onto the view area. we use a double-buffering system which means
			that we have a 'currently displayed' surface, and 'currently being drawn'
			surface. hence the 'swap' idea. in a single-buffering system we would see
			stuff being drawn one-after-the-other */
	while ( !glfwWindowShouldClose( g_window ) ) {
		/* wipe the drawing surface clear */
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		//glUseProgram(shader_programme);
		glBindVertexArray(vao);
		glLineWidth(1.5);
		/* draw points 0-3 from the currently bound VAO with current in-use shader */
		glDrawArrays(GL_POINTS, 0, totalCount);
		/* update other events like input handling */
		glfwPollEvents();
		/* put the stuff we've been drawing onto the display */
		glfwSwapBuffers( g_window );
	}
	
	/* close GL context and any other GLFW resources */
	glfwTerminate();
	return 0;
}
