/**************************************************************************************************\
*
* "universe.cpp" - Member function definitions for class universe ("universe.h")
*
*         Author - Tom McDonnell 2005
*
\**************************************************************************************************/

// Includes ////////////////////////////////////////////////////////////////////////////////////////

#include "universe.h"

#include <GL/glut.h>

#include <iostream>
#include <algorithm>

#include <assert.h>
#include <time.h>

// Global Variable Definitions /////////////////////////////////////////////////////////////////////

namespace bigBang
{
   extern rec2vector
   mainWinDim,     // Defined in main.cpp.
   halfMainWinDim; // Defined in main.cpp.
}

// Class Universe Member Function Definitions //////////////////////////////////////////////////////

namespace bigBang
{
   /*
    *
    */
   void universe::init(void)
   {
      unsigned int i, j;
      bool placed;

//      srand((unsigned)time(NULL));

      pVectorOld.resize(nParticles);
      pVectorNew.resize(nParticles);

      for (i = 0; i < pVectorOld.size(); ++i)
      {
         placed = false;

         while (!placed)
         {
            // Initialise particle positions to random within range [partRadius, max - partRadius].
            pVectorOld[i].pos.x = random(particleRadius, mainWinDim.x - particleRadius);
            pVectorOld[i].pos.y = random(particleRadius, mainWinDim.y - particleRadius);

            placed = true;

            // Detect collisions in initial placement.
            for (j = 0; j < i; ++j)
            {
               if (detectParticleOverlap(pVectorOld[i], pVectorOld[j]))
               {
                  placed = false;
               }
            }
         }

         // Initialise particle velocities to random within range [maxInitialV, -maxInitialV].
         pVectorOld[i].vel.x = random(-maxInitialV, maxInitialV);
         pVectorOld[i].vel.y = random(-maxInitialV, maxInitialV);
      }
/*
      // Override with test values.
      pVectorOld[0].pos = rec2vector(halfMainWinDim.x + 50.0,            100.0);
      pVectorOld[0].vel = rec2vector(                    0.0,              4.0);
      pVectorOld[1].pos = rec2vector(halfMainWinDim.x       , halfMainWinDim.y);
      pVectorOld[1].vel = rec2vector(                    0.0,             -2.0);
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

      for (unsigned int i = 0; i < nParticles; ++i)
      {
         glPushMatrix();
         glTranslatef(pVectorNew[i].pos.x, pVectorNew[i].pos.y, 0.0);
         glutWireSphere(particleRadius, 20, 10);
         glPopMatrix();
      }
   }

   /*
    *
    */
   void universe::simulate(double timeStep)
   {
      using std::cout;
      using std::endl;

      moveParticles(timeStep);
      bounceParticles();

      // Detect and resolve particle overlap.
      int p1No, p2No;
      double t = findTimeSince1stCollision(p1No, p2No);

      if (t > timeStep)
      {
         // No collisions have occurred.
         cout << "Collision occurred before timestep began!" << endl
              << "t        = " << t        << endl
              << "timeStep = " << timeStep << endl;
         exit(1);
      }

      if (t > 0.0)
      {
         // One or more collisions have occurred.

         //cout << "Pre -backtrack distance: "
         //     << distance(pVectorNew[p1No].pos, pVectorNew[p2No].pos)
         //     << endl;

         backTrack(-t); // Backtrack to time of 1st collision.

         //cout << "Post-backtrack distance: "
         //     << distance(pVectorNew[p1No].pos, pVectorNew[p2No].pos)
         //     << endl;

         // Calculate post-collision velocities for p1 & p2.
         collideParticles(pVectorNew[p1No], pVectorNew[p2No]);
         bounceParticles();
         pVectorOld.swap(pVectorNew);

         if (timeStep - t < 0.0)
         {
            //cout << " timeStep     = " << timeStep     << endl
            //     << " t            = " << t            << endl
            //     << " timeStep - t = " << timeStep - t << endl;
         }

         simulate(timeStep - t); // Simulate for remainder of timeStep.
      }
      else
      {
         // No collisions have occurred, or all have been dealt with.
         // New state of universe is legal so we can swap new and old particle arrays.
         pVectorOld.swap(pVectorNew);
      }
   }

   /*
    * Calculate new position and velocity for all particles,
    * resulting from their acceleration due to the gravitational
    * force between them and all other particles over the small period timeStep.
    * Data from the OLD array is used to produce new data for particle i in the NEW array.
    * This may result in overlap, as collisions are not detected here.
    */
   void universe::moveParticles(const double &timeStep)
   {
      assert(timeStep >= 0.0);

      if (timeStep < 0.0)
      {
         std::cout << "in moveParticles(), timeStep = " << timeStep << std::endl;
         exit(1);
      }

      for (unsigned int i = 0; i < nParticles; ++i)
      {
         rec2vector acc = rec2vector(0, 0);

// TODO: Figure out how to add touching particles to vector.
         std::vector<particle> particlesTouchingI;

         // Sum accelerations of particle i due to particles j.
         for (unsigned int j = 0; j < pVectorOld.size(); ++j)
         {
            if (i == j)
            {
               continue;
            }

            // Particles in pVectorOld are guaranteed not to be overlapping, but they may be very
            // close.  Remember the particles that are very close so that the component of the net
            // gravitational force on particle I in the direction of very close particles can be
            // made zero.  That way very frequent low speed collisions can be avoided.
            if (particlesAreTouching(pVectorOld[i], pVectorOld[j]))
            {
               particlesTouchingI[] = pVectorOld[j];
            }

            acc += convToRec(gravForce(pVectorOld[i], pVectorOld[j]) / particleMass);
         }

         // Set new velocity & position of particle i.
         pVectorNew[i].vel = pVectorOld[i].vel + acc * timeStep;
         pVectorNew[i].pos = pVectorOld[i].pos + pVectorNew[i].vel * timeStep;
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
      bool   collisionOccurred = false;
      double longestTime       = 0.0, t;

      for (unsigned int i = 0; i < pVectorNew.size(); ++i)
      {
         for (unsigned int j = 0; j < pVectorNew.size(); ++j)
         {
            if
            (
               i != j && detectParticleOverlap(pVectorNew[i], pVectorNew[j]) &&
               ((t = findTimeSinceCollision(i, j)) > longestTime)
            )
            {
               collisionOccurred = true;
               longestTime       =    t;
               p1No              =    i;
               p2No              =    j;
            }
         }
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
      using std::cout;
      using std::endl;

      particle &p1 = pVectorNew[p1No];
      particle &p2 = pVectorNew[p2No];

      assert(distance(p1.pos, p2.pos) < 2.0 * particleRadius); // Particles should be overlapping.

      double
      dvx = p2.vel.x - p1.vel.x        ,
      dvy = p2.vel.y - p1.vel.y        ,
      dx  = p2.pos.x - p1.pos.x        ,
      dy  = p2.pos.y - p1.pos.y        ,
      a   = pow(dvx, 2) + pow(dvy, 2)  ,
      b   = 2.0 * (dx * dvx + dy * dvy),
      c   = pow(dx, 2) + pow(dy, 2) - pow(2.0 * particleRadius, 2);

      double s1, s2;
      solvePolynomial(a, b, c, s1, s2);

      //cout << "s1, s2: " << s1 << ", " << s2 << endl;

      double d1 = distance(p1.pos                , p2.pos                ),
             d2 = distance(p1.pos + 1e-5 * p1.vel, p2.pos + 1e-5 * p2.vel);

      if (d2 > d1)
      {
         cout << "Particles should not be colliding."   << endl
              << "d1: "        << d1      <<               endl
              << "d2: "        << d2      <<               endl
              << "d2 - d1: "   << d2 - d1 <<               endl
              << "Solutions: " << s1      << ", " << s2 << endl;

         exit(1);
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
      using std::cout;
      using std::endl;

      assert(1.99 * particleRadius < distance(p1.pos, p2.pos) &&
                                     distance(p1.pos, p2.pos) < 2.01 * particleRadius);

      double
      PreSumVx = p1.vel.x + p2.vel.x,
      PreSumVy = p1.vel.y + p2.vel.y;

      // Create a unit vectors in direction p1->p2 and p2->p1.
      rec2vector p1p2 = p2.pos - p1.pos; p1p2 /= (2.0 * particleRadius);
      rec2vector p2p1 = -p1p2;

      // Calculate component of particles velocities in direction of impact.
      double
      p1Impact = vectDotProduct(p1.vel, p1p2),
      p2Impact = vectDotProduct(p2.vel, p2p1);

      // Calculate rebound velocity magnitude.
      double
      impactVMag  = p1Impact   + p2Impact,
      reboundVMag = impactVMag * reboundCoeff;

      // Update p1 and p2 velocities.
      p1.vel += reboundVMag * p2p1;
      p2.vel += reboundVMag * p1p2;

      double
      PostSumVx = p1.vel.x + p2.vel.x,
      PostSumVy = p1.vel.y + p2.vel.y;

      if (PreSumVx != PostSumVx || PreSumVy != PostSumVy)
      {
         //cout << "Sum vx Pre  = " << PreSumVx  << endl
         //     << "       Post = " << PostSumVx << endl
         //     << "Sum vy Pre  = " << PreSumVy  << endl
         //     << "       Post = " << PostSumVy << endl
         //     << endl;
      }
   }
}

/*******************************************END*OF*FILE*******************************************/
