/* pointcharge.c - Pointcharge at the origin */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

#define ffac (1/(4*gpt_pi*gpt_eps0))

struct pointcharge_info
{
  double mfac ;
} ;

static int pointcharge_sim(gptpar *par,double t,struct pointcharge_info *info) ;


void pointcharge_init(gptinit *init)
{
  struct pointcharge_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=1 )
    gpterror( "Syntax: %s(ECS,Q)\n", gptgetname(init) ) ;

  info = (struct pointcharge_info *)gptmalloc( sizeof(struct pointcharge_info) ) ;

  info->mfac = ffac*gptgetargdouble(init,1) ;

  gptaddEBelement( init, pointcharge_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}


static int pointcharge_sim(gptpar *par,double t,struct pointcharge_info *info)
{
  double r2, r3 ;

  r2=X*X+Y*Y+Z*Z ;
  if( r2==0.0 ) return( 0 ) ;
  r3=r2*sqrt(r2) ;

  EX=info->mfac*X/r3 ;
  EY=info->mfac*Y/r3 ;
  EZ=info->mfac*Z/r3 ;

  return( 1 ) ;
}
