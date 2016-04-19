/**************************************************************************************************\
*
* "main.cpp" - Main file of the 'Big Bang' program.
*              Contains the main() function plus initialisation and event handling functions.
*
*     Author - Tom McDonnell 2004
*
\**************************************************************************************************/

// Includes ////////////////////////////////////////////////////////////////////////////////////////

#include "lib_tom_cpp/vector.h"
#include "universe.h"

#include <GL/glut.h>

#include <iostream>

// Global Variable Definitions /////////////////////////////////////////////////////////////////////

namespace bigBang
{
   rec2vector
   mainWinDim(600, 600),
   halfMainWinDim(mainWinDim.x / 2, mainWinDim.y / 2);

   universe U;
}

// File Scope Function Declarations ////////////////////////////////////////////////////////////////

namespace
{
   void display(void);
   void reshape(int, int);
   void idle(void);
   void setupView(void);
}

// Main Function Definition ////////////////////////////////////////////////////////////////////////

/*
 *
 */
int main(int argc, char *argv[])
{
   using namespace bigBang;

   // Initialise glut.
   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
   
   // Create main window.
   glutInitWindowSize(mainWinDim.x, mainWinDim.y);
   glutInitWindowPosition(0, 0);
   glutCreateWindow("Big Bang Simulator");

   // Register functions for main window.
   glutDisplayFunc(display);
   glutReshapeFunc(reshape);
   glutIdleFunc(idle);

   U.init();

   glutMainLoop();

   return 0;
}

// File Scope Function Definitions /////////////////////////////////////////////////////////////////

namespace
{
   using namespace bigBang;

   /*
    *
    */
   void display(void)
   {
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      U.draw();

      setupView();

      glutSwapBuffers();
   }

   /*
    *
    */
   void reshape(int w, int h)
   {
      using namespace bigBang;

      mainWinDim.x   = w;
      mainWinDim.y   = h;
      halfMainWinDim = mainWinDim / 2; // keep as integer
      glutPostRedisplay();
   }

   /*
    *
    */
   void idle(void)
   {
      U.simulate(1.0);
     
      glutPostRedisplay();
   }

   /*
    * Set up modelview matrix for viewing x-y plane from -tve z direction.
    */
   void setupView(void)
   {
      using namespace bigBang;

      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      glOrtho
      (
         -halfMainWinDim.x,  halfMainWinDim.x, // Left  , right.
         -halfMainWinDim.y,  halfMainWinDim.y, // Bottom, top  .
         0.0              ,  200.0             // Near  , far  .
      );

      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
      gluLookAt
      (
         halfMainWinDim.x, halfMainWinDim.y, 100.0, // Look from.
         halfMainWinDim.x, halfMainWinDim.y,   0.0, // Look at  .
         0.0             ,              1.0,   0.0  // Up vector.
      );
   }
}

/*******************************************END*OF*FILE********************************************/
