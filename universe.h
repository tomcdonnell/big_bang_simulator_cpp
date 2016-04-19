/**************************************************************************************************\
*
* "universe.h" -
*
*       Author - Tom McDonnell 2005
*
\**************************************************************************************************/

#ifndef UNIVERSE_H
#define UNIVERSE_H

// Includes ////////////////////////////////////////////////////////////////////////////////////////

#include "lib_tom_cpp/vector.h"

#include <fstream>
#include <vector>

// Class Definitions ///////////////////////////////////////////////////////////////////////////////

namespace bigBang
{
   using namespace TomsLibVector;

   class particle
   {
    public:
      rec2vector pos, vel;
   };

   /*
    *
    */
   class universe
   {
    public:

      /*
       * Constructor.
       */
      universe
      (
         const double _G                 =  1.0e3,
         const int    _nParticles        =  3    ,
         const double _particleRadius    = 50.0  ,
         const double _particleMass      =  1.0  , 
         const double _maxInitialV       =  0.0  ,
         const double _reboundCoeff      =  0.9  ,
         const double _touchingTolerance =  1.0
      )
      :
      nParticles(_nParticles)        ,
      G(_G)                          ,
      particleRadius(_particleRadius),
      particleMass(_particleMass)    ,
      maxInitialV(_maxInitialV)      ,
      reboundCoeff(_reboundCoeff)    ,
      touchingTolerance(_touchingTolerance)
      {}

      void init(void);
      void draw(void);
      void simulate(double timeStep);

    private:

      inline pol2vector             gravForce(const particle &, const particle &);
      inline bool       detectParticleOverlap(const particle &, const particle &);
      inline bool        particlesAreTouching(const particle &, const particle &);

      inline void moveParticles(const double &timeStep);
      inline void     backTrack(const double &timeStep);

      inline void bounceParticles(void);
      inline void   wrapParticles(void);

      double findTimeSince1stCollision(int &p1No, int &p2No);
      double findTimeSinceCollision(const int &p1No, const int &p2No);

      void collideParticles(particle &, particle &);

      const unsigned int nParticles;

      const double
      G,                 // (Nm^2 / kg^2) universal gravitational constant.
                         //   Chosen arbitrarily, but for consistent velocities,
                         //   needs to be inv-proportional to nParticles
                         //   (fewer particles need stronger gravity).
      particleRadius,    // (m) (1 pixel = 1m)
      particleMass  ,    // (kg)
      maxInitialV   ,    // (m/s)
      reboundCoeff  ,    // ratio: (post-collision energy) / (pre-collision energy)
      touchingTolerance; // (m) (1 pixel = 1m).

      std::vector<particle>
      pVectorOld, // Old array is always in consistent state (ie. no overlap).
      pVectorNew; // New array is used during simulation to ensure that all inconsistencies are
                  // removed before old becomes new and vice-versa.
   };
}

#endif

// Inline Member Function Definitions /////////////////////////////////////////////////////////////

namespace bigBang
{
   using namespace TomsLibVector;

   /*
    * Calculate the gravitational force on p1 due to the
    * gravitational field interaction between p1 and p2.
    */
   inline pol2vector universe::gravForce(const particle &p1, const particle &p2)
   {
      pol2vector f = convToPol(p2.pos - p1.pos);

      f.setR((G * 2.0 * particleMass) / pow(f.getR(), 2.0));

      return f;
   }

   /*
    *
    */
   inline void universe::backTrack(const double &timeStep)
   {
      assert(timeStep <= 0.0);

      for (unsigned int i = 0; i < nParticles; ++i)
      {
         pVectorNew[i].pos = pVectorNew[i].pos + pVectorNew[i].vel * timeStep;
      }
   }

   extern rec2vector mainWinDim; // defined in "main.cpp"

   /*
    * If particles positions are outside screen boundaries,
    * Calculate new positions mirrored about the border line(s) that was(were) exceeded.
    * (result is particles bounce off walls)
    */
   inline void universe::bounceParticles(void)
   {
      for (unsigned int i = 0; i < nParticles; ++i)
      {
         const double &r = particleRadius;
         rec2vector   &p = pVectorNew[i].pos,
                      &v = pVectorNew[i].vel;

         // Bounce x coordinate and velocity.
              if (p.x > mainWinDim.x - r) {p.x = 2.0 * (mainWinDim.x - r) - p.x; v.x = -v.x;}
         else if (p.x <                r) {p.x = 2.0 *                 r  - p.x; v.x = -v.x;}

         // Bounce y coordinate and velocity.
              if (p.y > mainWinDim.y - r) {p.y = 2.0 * (mainWinDim.y - r) - p.y; v.y = -v.y;}
         else if (p.y <                r) {p.y = 2.0 *                 r  - p.y; v.y = -v.y;}
      }
   }

   /*
    * Alternative to bounce().
    * If particles positions are outside screen boundaries,
    * Calculate new positions by wrapping coordinates through opposite walls.
    * (result is particles disappear through one wall, and reappear at the opposite one)
    */
   inline void universe::wrapParticles(void)
   {
      for (unsigned int i = 0; i < nParticles; ++i)
      {
         rec2vector   &p = pVectorNew[i].pos;

         // Wrap x coordinate.
              if (p.x > mainWinDim.x) p.x -= mainWinDim.x;
         else if (p.x <          0.0) p.x += mainWinDim.x;

         // Wrap y coordinate.
              if (p.y > mainWinDim.y) p.y -= mainWinDim.y;
         else if (p.y <          0.0) p.y += mainWinDim.y;
      }
   }

   /*
    * Detect overlap between particles p1 & p2.
    */
   inline bool universe::detectParticleOverlap(const particle &p1, const particle &p2)
   {   
      const double &r = particleRadius;

      return
      (
         p1.pos.x + r > p2.pos.x - r &&
         p1.pos.x - r < p2.pos.x + r &&
         p1.pos.y + r > p2.pos.y - r &&
         p1.pos.y - r < p2.pos.y + r &&
         distance(p1.pos, p2.pos) < 2.0 * particleRadius
      );
   }

   /*
    * If the distance between p1 & p2 is less than the touching tolerance, return true else false.
    */
   inline bool universe::particlesAreTouching(const particle &p1, const particle &p2)
   {   
      const double &r = particleRadius + touchingTolerance;

      return
      (
         p1.pos.x + r > p2.pos.x - r &&
         p1.pos.x - r < p2.pos.x + r &&
         p1.pos.y + r > p2.pos.y - r &&
         p1.pos.y - r < p2.pos.y + r &&
         distance(p1.pos, p2.pos) < 2.0 * particleRadius
      );
   }
}

/*******************************************END*OF*FILE********************************************/
