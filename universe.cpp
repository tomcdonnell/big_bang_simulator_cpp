/*************************************************************************************************\
*                                                                                                 *
*  "universe.cpp" - Member function definitions for class universe ("universe.h")                 *
*                                                                                                 *
*          Author - Tom McDonnell 2005                                                            *
*                                                                                                 *
\*************************************************************************************************/

// INCLUDES ///////////////////////////////////////////////////////////////////////////////////////

#include "universe.h"

#include <GL/glut.h>

#include <iostream>
#include <algorithm>

#include <assert.h>
#include <time.h>

// GLOBAL VARIABLE DEFINITIONS ////////////////////////////////////////////////////////////////////

namespace bigBang
{

 extern rec2vector mainWinDim,     // defined in main.cpp
                   halfMainWinDim; // defined in main.cpp

}

// CLASS UNIVERSE MEMBER FUNCTION DEFINITIONS /////////////////////////////////////////////////////

namespace bigBang
{

 /*
  *
  */
 void universe::init(void)
 {
    int i, j;
    bool placed;

    //srand((unsigned)time(NULL));

    pVectorOld.resize(nParticles);
    pVectorNew.resize(nParticles);

    for (i = 0; i < pVectorOld.size(); ++i)
    {
       placed = false;
       while (!placed)
       {
          // initialise particle positions to random within range [partRadius, max - partRadius]
          pVectorOld[i].pos.x = random(particleRadius, mainWinDim.x - particleRadius);
          pVectorOld[i].pos.y = random(particleRadius, mainWinDim.y - particleRadius);

          placed = true;

          // detect collisions in initial placement
          for (j = 0; j < i; ++j)
            if (detectParticleOverlap(pVectorOld[i], pVectorOld[j]))
              placed = false;
       }

       // initialise particle velocities to random within range [maxInitialV, -maxInitialV]
       pVectorOld[i].vel.x = random(-maxInitialV, maxInitialV);
       pVectorOld[i].vel.y = random(-maxInitialV, maxInitialV);
    }
/*
    // override with test values
    pVectorOld[0].pos = rec2vector(halfMainWinDim.x + 50.0,            100.0);
    pVectorOld[0].vel = rec2vector(             0.0,              4.0);
    pVectorOld[1].pos = rec2vector(halfMainWinDim.x, halfMainWinDim.y);
    pVectorOld[1].vel = rec2vector(             0.0,             -2.0);
*/
    glutPostRedisplay();
 }

 /*
  * Draw universe in new position.
  */
 void universe::draw(void)
 {
    glColor3f(1.0, 1.0, 1.0);

    glMatrixMode(GL_MODELVIEW);

    //glBegin(GL_POINTS);
      for (int i = 0; i < nParticles; ++i)
      {
         //glVertex3f(pArrayOld[i].pos.x, pArrayOld[i].pos.y, 0.0);
         glPushMatrix();
         glTranslatef(pVectorNew[i].pos.x, pVectorNew[i].pos.y, 0.0);
         glutWireSphere(particleRadius, 20, 10);
         glPopMatrix();
      }
   //glEnd();
 }

 /*
  *
  */
 void universe::simulate(double timeStep)
 {
    moveParticles(timeStep);
    bounceParticles();     // (particles may have travelled off screen)

    // detect and resolve particle overlap
    int p1No, p2No;
    double t = findTimeSince1stCollision(p1No, p2No);

    if (t > timeStep) cout << "t        = " << t        << endl
                           << "timeStep = " << timeStep << endl;

    if (t > 0.0)
    {
       // one or more collisions have occurred
/*
       cout << "Pre -backtrack distance: "
            << distance(pVectorNew[p1No].pos, pVectorNew[p2No].pos)
            << endl;
*/
       backTrack(-t);                                    // backtrack to time of 1st collision
/*
       cout << "Post-backtrack distance: "
            << distance(pVectorNew[p1No].pos, pVectorNew[p2No].pos)
            << endl;
*/
       collideParticles(pVectorNew[p1No], pVectorNew[p2No]); // calculate post-collision velocities
                                                             //   for p1 & p2
       bounceParticles();
       pVectorOld.swap(pVectorNew);

       if ((timeStep - t) < 0.0)
         cout << " timeStep     = " << timeStep     << endl
              << " t            = " << t            << endl
              << " timeStep - t = " << timeStep - t << endl;

       simulate(timeStep - t);                               // simulate for remainder of timeStep
    }
    else
    {
       // no collisions have occurred, or all have been dealt with

       // new state of universe is legal so we
       //   can swap new and old particle arrays
       pVectorOld.swap(pVectorNew);
    }
 }

 /*
  * Check every pair of particles for possible collision
  * while keeping track of: longest time since collision seen so far, and
  *                         particle numbers involved in that collision.
  * The collision we are left with is the one sought.
  * Returns the time since the 1st collision, and leaves the particle numbers in p1No & p2No.
  */
 double universe::findTimeSince1stCollision(int &p1No, int &p2No)
 {
    bool collisionOccurred = false;

    double longestTime = 0.0, t;
    for (int i = 0; i < pVectorNew.size(); ++i)
      for (int j = 0; j < pVectorNew.size(); ++j)
        if (i != j && detectParticleOverlap(pVectorNew[i], pVectorNew[j]))
          if ((t = findTimeSinceCollision(i, j)) > longestTime)
          {
             collisionOccurred = true;
             longestTime       =    t;
             p1No              =    i;
             p2No              =    j;
          }

    if (collisionOccurred) return longestTime;
    else                   return        -1.0;
 }

 /*
  * Assuming p1 has been found to be overlapping p2, find points
  * where p1 and p2 must have been when initial contact took place.
  * Return the time elapsed since the collision took place.
  */
 double universe::findTimeSinceCollision(const int &p1No, const int &p2No)
 {
    particle &p1 = pVectorNew[p1No];
    particle &p2 = pVectorNew[p2No];

    assert(distance(p1.pos, p2.pos) < 2.0 * particleRadius); // particles should be overlapping

    double dvx = p2.vel.x - p1.vel.x,
           dvy = p2.vel.y - p1.vel.y,
           dx  = p2.pos.x - p1.pos.x,
           dy  = p2.pos.y - p1.pos.y,
           a = pow(dvx, 2) + pow(dvy, 2),
           b = 2.0 * (dx * dvx + dy * dvy),
           c = pow(dx, 2) + pow(dy, 2) - pow(2.0 * particleRadius, 2);

    double s1, s2;
    solvePolynomial(a, b, c, s1, s2);

    double d1 = distance(p1.pos                , p2.pos                ),
           d2 = distance(p1.pos + 1e-5 * p1.vel, p2.pos + 1e-5 * p2.vel);
    if (d2 > d1)
    {
       cout << "Particles should not be colliding. d2 - d1 = " << d2 - d1 << endl
            << "Solutions: " << s1 << ", " << s2 << endl;

       //return 0.0;
    }

         if (s1 < 0.0 && s2 > 0.0) return -s1;
    else if (s2 < 0.0 && s1 > 0.0) return -s2;
    else
    {
       cout << "Solutions: " << s1 << ", " << s2 << endl;
       error("Unexpected result in findInitContactPosition().");
       return 0.0;
    }
 }

 /*
  * Assuming p1 & p2 are have just collided,
  * (p1 & p2 are touching, and their according to their velocities, are moving towards each other)
  * update the velocities of p1 & p2 to their post-collision values
  * (they will then be moving away from each other).
  */
 void universe::collideParticles(particle &p1, particle &p2)
 {
    assert(   1.99 * particleRadius < distance(p1.pos, p2.pos)
           &&                         distance(p1.pos, p2.pos) < 2.01 * particleRadius);

    double PreSumVx = p1.vel.x + p2.vel.x,
           PreSumVy = p1.vel.y + p2.vel.y;

    // create a unit vectors in direction p1->p2 and p2->p1
    rec2vector p1p2 = p2.pos - p1.pos; p1p2 /= (2.0 * particleRadius);
    rec2vector p2p1 = -p1p2;

    // calculate component of particles velocities in direction of impact
    double p1Impact = vectDotProduct(p1.vel, p1p2),
           p2Impact = vectDotProduct(p2.vel, p2p1);

    // calculate rebound velocity magnitude
    double impactVMag  = p1Impact + p2Impact,
           reboundVMag = impactVMag * reboundCoeff;

    // Update p1 and p2 velocities
    p1.vel += reboundVMag * p2p1;
    p2.vel += reboundVMag * p1p2;

    double PostSumVx = p1.vel.x + p2.vel.x,
           PostSumVy = p1.vel.y + p2.vel.y;

    if (PreSumVx != PostSumVx || PreSumVy != PostSumVy)
      cout << "Sum vx Pre  = " << PreSumVx  << endl
           << "       Post = " << PostSumVx << endl
           << "Sum vy Pre  = " << PreSumVy  << endl
           << "       Post = " << PostSumVy << endl
           << endl;
 }

} // end namespace bigBang

/*******************************************END*OF*FILE*******************************************/
