/* rmax.c - Kill particle when R is too large */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "elem.h"

struct rmax_info
{
  double rmax2 ;
  double zlen ;
} ;

static int rmax_sim(gptpar *par,double t,struct rmax_info *info) ;


void rmax_init(gptinit *init)
{
  struct rmax_info *info ;
  int numarg ;

  gptbuildECS( init ) ;

  numarg = gptgetargnum(init) ;

  if( numarg!=1 && numarg!=2 )
    gpterror( "Syntax: %s(ECS,R,[L])\n", gptgetname(init) ) ;

  info = (struct rmax_info *)gptmalloc( sizeof(struct rmax_info) ) ;

  info->rmax2 = SQR(gptgetargdouble(init,1)) ;
  if( numarg==2 )
    info->zlen = gptgetargdouble(init,2)/2 ;
  else
    info->zlen = DBL_MAX ;

  gptaddEBelement( init, rmax_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}


static int rmax_sim(gptpar *par,double t,struct rmax_info *info)
{
  if( fabs(Z) > info->zlen ) return( 0 ) ;
  if( SQR(X)+SQR(Y) <= info->rmax2 ) return( 1 ) ;

  gptremoveparticle(par) ;
  return( 1 ) ;
}
