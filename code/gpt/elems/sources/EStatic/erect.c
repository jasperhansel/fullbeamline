/* erect.c - Rectangular local E-field in the z-direction */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

struct erect_info
{
  double a,b,L ;
  double E ;
} ;

static int erect_sim(gptpar *par,double t,struct erect_info *info) ;


void erect_init(gptinit *init)
{
  struct erect_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=4 )
    gpterror( "Syntax: %s(ECS,a,b,L,E)\n", gptgetname(init) ) ;

  info = (struct erect_info *)gptmalloc( sizeof(struct erect_info) ) ;

  info->a = gptgetargdouble(init,1)/2 ;
  info->b = gptgetargdouble(init,2)/2 ;
  info->L = gptgetargdouble(init,3)/2 ;
  info->E = gptgetargdouble(init,4) ;

  gptaddEBelement( init, erect_sim, gptfree, GPTELEM_LOCAL, info ) ;
}


static int erect_sim(gptpar *par,double t,struct erect_info *info)
{
  if( fabs(X)>info->a ||
      fabs(Y)>info->b ||
      fabs(Z)>info->L ) return( 0 ) ;

  EZ = info->E ;

  return( 1 ) ;
}
