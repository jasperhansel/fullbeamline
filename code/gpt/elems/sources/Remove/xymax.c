/* xymax.c - Kill particle when |x| or |y| coordinate is too large */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "elem.h"

struct xymax_info
{
  double a ;
  double b ;
  double zlen ;
} ;

static int xymax_sim(gptpar *par,double t,struct xymax_info *info) ;


void xymax_init(gptinit *init)
{
  struct xymax_info *info ;
  int numarg ;

  gptbuildECS( init ) ;

  numarg = gptgetargnum(init) ;

  if( numarg!=2 && numarg!=3 )
    gpterror( "Syntax: %s(ECS,a,b,[L])\n", gptgetname(init) ) ;

  info = (struct xymax_info *)gptmalloc( sizeof(struct xymax_info) ) ;

  info->a = gptgetargdouble(init,1)/2 ;
  info->b = gptgetargdouble(init,2)/2 ;
  if( numarg==3 )
    info->zlen = gptgetargdouble(init,3)/2 ;
  else
    info->zlen = DBL_MAX ;

  gptaddEBelement( init, xymax_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}


static int xymax_sim(gptpar *par,double t,struct xymax_info *info)
{
  if( fabs(Z) > info->zlen ) return( 0 ) ;
  if( fabs(X) <= info->a && fabs(Y) <= info->b ) return( 1 ) ;

  gptremoveparticle(par) ;
  return( 1 ) ;
}
