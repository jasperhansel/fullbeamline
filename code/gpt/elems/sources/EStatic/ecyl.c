/* ecyl.c - Cylindric local E-field in the z-direction */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

struct ecyl_info
{
  double field ;
  double zlen ;
  double radius2 ;
} ;

static int ecyl_sim(gptpar *par,double t,struct ecyl_info *info) ;


void ecyl_init(gptinit *init)
{
  struct ecyl_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=3 )
    gpterror( "Syntax: %s(ECS,R,L,E)\n", gptgetname(init) ) ;

  info = (struct ecyl_info *)gptmalloc( sizeof(struct ecyl_info) ) ;

  info->radius2 = SQR(gptgetargdouble(init,1)) ;
  info->zlen    = gptgetargdouble(init,2)/2 ;
  info->field   = gptgetargdouble(init,3) ;

  gptaddEBelement( init, ecyl_sim, gptfree, GPTELEM_LOCAL, info ) ;
}


static int ecyl_sim(gptpar *par,double t,struct ecyl_info *info)
{
  if( fabs( Z ) > info->zlen ) return( 0 ) ;
  if( SQR(X)+SQR(Y) > info->radius2 ) return( 0 ) ;

  EZ = info->field ;

  return( 1 ) ;
}
