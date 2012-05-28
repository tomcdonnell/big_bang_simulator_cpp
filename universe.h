/*************************************************************************************************\
*                                                                                                 *
*  "universe.h" -                                                                                 *
*                                                                                                 *
*        Author - Tom McDonnell 2005                                                              *
*                                                                                                 *
\*************************************************************************************************/

#ifndef UNIVERSE_H
#define UNIVERSE_H

// INCLUDES ///////////////////////////////////////////////////////////////////////////////////////

#include "lib_tom_cpp/vector.h"

#include <fstream>
#include <vector>

// CLASS DEFINITIONS //////////////////////////////////////////////////////////////////////////////

namespace bigBang
{
 using namespace TomsLibVector;

 class particle
 {
  public:
   
    rec2vector pos,
               vel;
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
    universe(const double _G              =  1.0e3, const int    _nParticles     = 10 ,
             const double _particleRadius = 50.0  , const double _particleMass   = 1.0, 
             const double _maxInitialV    =  4.0  , const double _reboundCoeff   = 1.0 )
    : nParticles(_nParticles)        , G(_G)                      ,
      particleRadius(_particleRadius), particleMass(_particleMass),
      maxInitialV(_maxInitialV)      , reboundCoeff(_reboundCoeff)
    {}

    void init(void);
    void draw(void);
    void simulate(double timeStep);

    //void writeParticlePositionsToFile(ofstream);

  private:

    inline pol2vector             gravForce(const particle &, const particle &);
    inline bool       detectParticleOverlap(const particle &, const particle &);

    inline void moveParticles(const double &timeStep);
    inline void     backTrack(const double &timeStep);

    inline void bounceParticles(void);
    inline void   wrapParticles(void);

    double findTimeSince1stCollision(int &p1No, int &p2No);
    double findTimeSinceCollision(const int &p1No, const int &p2No);

    void collideParticles(particle &, particle &);

    const unsigned int nParticles;

    const double G,              // (Nm^2 / kg^2) universal gravitational constant.
                                 //   Chosen arbitrarily, but for consistent velocities,
                                 //   needs to be inv-proportional to nParticles
                                 //   (fewer particles need stronger gravity).
                 particleRadius, // (m) (1 pixel = 1m)
                 particleMass,   // (kg)
                 maxInitialV,    // (m/s)
                 reboundCoeff;   // ratio: (post-collision energy) / (pre-collision energy)

    std::vector<particle> pVectorOld, // old array is always in consistent state (ie. no overlap).
                          pVectorNew; // new array is used during simulation to ensure that all
                                      //   inconsistencies are removed before old becomes new
                                      //   and vice-versa.
   };

} // end namespace bigBang

#endif

// INLINE MEMBER FUNCTION DEFINITIONS /////////////////////////////////////////////////////////////

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

    //if (f.getR() < particleRadius)
      // zero force if particles are overlapping
      f.setR(0.0);
    //else
      //f.setR((G * 2.0 * particleMass) / pow(f.getR(), 2.0));

    return f;
 }

 /*
  * Calculate new position and velocity for all particles,
  * resulting from their acceleration due to the gravitational
  * force between them and all other particles over the small period timeStep.
  * Data from the OLD array is used to produce new data for particle i in the NEW array.
  * This may result in overlap, as collisions are not detected here.
  */
 inline void universe::moveParticles(const double &timeStep)
 {
    //assert(timeStep >= 0.0);
    if (timeStep < 0.0)
    {
       std::cout << "in moveParticles(), timeStep = " << timeStep << std::endl;
       exit(1);
    }

    for (unsigned int i = 0; i < nParticles; ++i)
    {
       rec2vector acc = rec2vector(0, 0);

       // sum accelerations of particle i due to particles j
       for (unsigned int j = 0; j < pVectorOld.size(); ++j)
         if (i != j)
           acc += convToRec(gravForce(pVectorOld[i], pVectorOld[j]) / particleMass);

       // set new velocity & position of particle i
       pVectorNew[i].vel = pVectorOld[i].vel + acc * timeStep;
       pVectorNew[i].pos = pVectorOld[i].pos + pVectorNew[i].vel * timeStep;
    }
 }

 /*
  *
  */
 inline void universe::backTrack(const double &timeStep)
 {
    assert(timeStep <= 0.0);

    for (unsigned int i = 0; i < nParticles; ++i)
      pVectorNew[i].pos = pVectorNew[i].pos + pVectorNew[i].vel * timeStep;
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

       // bounce x coordinate and velocity
            if (p.x > mainWinDim.x - r) {p.x = 2.0 * (mainWinDim.x - r) - p.x; v.x = -v.x;}
       else if (p.x <                r) {p.x = 2.0 *                 r  - p.x; v.x = -v.x;}

       // bounce y coordinate and velocity
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

       // wrap x coordinate
            if (p.x > mainWinDim.x) p.x -= mainWinDim.x;
       else if (p.x <          0.0) p.x += mainWinDim.x;

       // wrap y coordinate
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

    if (   p1.pos.x + r > p2.pos.x - r
        && p1.pos.x - r < p2.pos.x + r
        && p1.pos.y + r > p2.pos.y - r
        && p1.pos.y - r < p2.pos.y + r
        && distance(p1.pos, p2.pos) < 2.0 * particleRadius)
      return true;
    else
      return false;
 }

} // end namespace bigBang

/*******************************************END*OF*FILE*******************************************/
