/* magdipole.c: Magnetic mathematical dipole */

#include <stdio.h>
#include <math.h>
#include "elem.h"

/* Info structure containing all relevant parameters for this element */
struct magdipole_info
{
  double m ;
} ;

/* Forward declaration of the routine calculating the electromagnetic fields */
static int magdipole_sim(gptpar *par,double t,struct magdipole_info *info) ;

/* Initialization routine */
void magdipole_init(gptinit *init)
{
  struct magdipole_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=1 )
    gpterror( "Syntax: %s(ECS,m)\n", gptgetname(init) ) ;

  info = (struct magdipole_info *)gptmalloc( sizeof(struct magdipole_info) ) ;
  info->m        = gptgetargdouble(init,1) ;

  gptaddEBelement( init, magdipole_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}


/* The following routine calculates the electromagnetic fields */
static int magdipole_sim(gptpar *par,double t,struct magdipole_info *info)
{
  double m  = info->m ;
  double r2 = X*X+Y*Y+Z*Z ;
  double r  = sqrt(r2) ;
  double fac ;

  if( r==0 )
  {
    gptwarning( "Infinite field in magdipole") ;
    return( 0 ) ;
  }

  fac= (gpt_mu0*m)/(4*gpt_pi*r*r2*r2) ;

  BX = 3*fac*Z*X ;
  BY = 3*fac*Z*Y ;
  BZ = fac*(3*Z*Z-r2) ;

  return( 1 ) ;
}
