//g++ -o version.exe version.cpp libglfw3dll.a libglew32.dll.a –I include –lglew32 -lglfw3 –lgdi32 –lopengl32

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <stdio.h>

int main (){
	if (!glfwInit ()){
		fprintf (stderr, "ERROR: could not start GLFW3\n");
		return 1;
	}
	GLFWwindow* window = glfwCreateWindow( 640, 480, "Hello Triangle", NULL, NULL );
	if ( !window ) {
		fprintf( stderr, "ERROR: could not open window with GLFW3\n" );
		glfwTerminate();
		return 1;
	}
	glfwMakeContextCurrent( window );
	/* start GLEW extension handler */
	glewExperimental = GL_TRUE;
	glewInit();

	/* get version info */
	const GLubyte* renderer = glGetString( GL_RENDERER ); /* get renderer string */
	const GLubyte* version = glGetString( GL_VERSION );	 /* version as a string */
	printf( "Renderer: %s\n", renderer );
	printf( "OpenGL version supported %s\n", version );
	
	/* tell GL to only draw onto a pixel if the shape is closer to the viewer
	than anything already drawn at that pixel */
	glEnable( GL_DEPTH_TEST ); /* enable depth-testing */
	/* with LESS depth-testing interprets a smaller depth value as meaning "closer" */
	glDepthFunc( GL_LESS );
	
	/* close GL context and any other GLFW resources */
	glfwTerminate();
	return 0;
}