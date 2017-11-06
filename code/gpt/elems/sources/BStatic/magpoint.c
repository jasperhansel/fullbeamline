/* magpoint.c - Magnetic pointcharge at the origin */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

#define ffac (gpt_mu0/(4*gpt_pi))

struct magpoint_info
{
  double mfac ;
} ;

static int magpoint_sim(gptpar *par,double t,struct magpoint_info *info) ;


void magpoint_init(gptinit *init)
{
  struct magpoint_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=1 )
    gpterror( "Syntax: %s(ECS,Q)\n", gptgetname(init) ) ;

  info = (struct magpoint_info *)gptmalloc( sizeof(struct magpoint_info) ) ;

  info->mfac = ffac*gptgetargdouble(init,1) ;

  gptaddEBelement( init, magpoint_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}


static int magpoint_sim(gptpar *par,double t,struct magpoint_info *info)
{
  double r2, r3 ;

  r2=X*X+Y*Y+Z*Z ;
  if( r2==0.0 ) return( 0 );

  r3=r2*sqrt(r2) ;

  BX=info->mfac*X/r3 ;
  BY=info->mfac*Y/r3 ;
  BZ=info->mfac*Z/r3 ;

  return( 1 ) ;
}
