/* ehole.c: Uniform fields separated by a conducting plate with circular hole */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 *
 * Equations from 3d edition Jackson pp. 130.
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

struct ehole_info
{
  double a ;
  double E0 ;
  double E1 ;
} ;

static int ehole_sim(gptpar *par,double t,struct ehole_info *info) ;

/* Initialization routine */
void ehole_init(gptinit *init)
{
  struct ehole_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=3 )
    gpterror( "Syntax: %s(ECS,R,E1,E0)\n", gptgetname(init) ) ;

  info = (struct ehole_info *)gptmalloc( sizeof(struct ehole_info) ) ;

  info->a        = gptgetargdouble(init,1) ;
  info->E1       = gptgetargdouble(init,2) ;
  info->E0       = gptgetargdouble(init,3) ;

  gptaddEBelement( init, ehole_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}


/* The following routine calculates the electromagnetic fields */
static int ehole_sim(gptpar *par,double t,struct ehole_info *info)
{
  /* Copy of parameters in info structure for convenience */
  double a, E0, E1 ;
  double R, lambda ;
  double sm, sp ;
  double signz ;
  double mfac ;
  double Ez, Eror ;

  /* Retrieve parameters from info structure */
  a      = info->a ;
  E0     = info->E0 ;
  E1     = info->E1 ;

  /* Only calculate fringe fields when a!=0 */
  if( a>0 )
  {
    lambda = (Z*Z+X*X+Y*Y - a*a)/(a*a) ;
    R      = sqrt(lambda*lambda+4*Z*Z/(a*a)) ;

    sm = sqrt((R-lambda)/2) ;
    sp = sqrt((R+lambda)/2) ;

    signz = 0 ;
    if( Z<0 ) signz = -1 ;
    if( Z>0 ) signz = +1 ;

    Ez = 0 ;
    if( sm!=0.0 )
      Ez   += Z*(1/sm-sm) ;

    Eror = -sm ;
    if( sp!=0.0 )
    {
      Ez   += Z*fabs(Z)/(a*sp) ;
      Ez   -= a*R*atan(1/sp)*signz ;
      Eror += (fabs(Z)/a) * (1/(1/sp+sp)) ;
    }
    else
      Ez   -= a*R*gpt_pi/2 ;

    mfac = (E0-E1)/(a*gpt_pi*R) ;

    EZ = mfac*Ez ;
    EX = mfac*X*Eror ;
    EY = mfac*Y*Eror ;
  }

  if( Z<0 ) EZ += E1 ;
  if( Z>0 ) EZ += E0 ;

  /* Particle is always inside element */
  return( 1 ) ;
}
