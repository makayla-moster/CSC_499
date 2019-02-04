
// g++ -w -o monkm.exe CImg_Test.cpp -O2 -lgdi32

#include "CImg.h"
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

using namespace cimg_library;
using namespace std;


int main() {
  CImg<unsigned char> img(640,400,1,3);  // Define a 640x400 color image with 8 bits per color component.
  img.fill(0);                           // Set pixel values to 0 (color : black)
  unsigned char purple[] = { 255,0,255 };        // Define a purple color
  img.draw_text(100,100,"Hello World",purple); // Draw a purple "Hello world" at coordinates (100,100).
  img.display("My first CImg code");             // Display the image in a display window.
  return 0;
}
