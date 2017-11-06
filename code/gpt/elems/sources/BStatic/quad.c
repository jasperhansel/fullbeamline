/* quad.c - Quadrupole lens */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

struct quad_info
{
  double G ;
  double zlen ;
} ;

static int quad_sim(gptpar *par, double t, struct quad_info *info) ;


void quadrupole_init(gptinit *init)
{
  struct quad_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=2 )
    gpterror( "Syntax: %s(ECS,L,G)\n", gptgetname(init) ) ;

  info = (struct quad_info *)gptmalloc( sizeof(struct quad_info) ) ;

  info->zlen = gptgetargdouble(init,1)/2 ;
  info->G    = gptgetargdouble(init,2) ;

  gptaddEBelement( init, quad_sim, gptfree, GPTELEM_LOCAL, info ) ;
}


static int quad_sim(gptpar *par, double t, struct quad_info *info)
{
  if( fabs(Z) > info->zlen ) return( 0 ) ;

  BX = info->G * Y ;
  BY = info->G * X ;

  return( 1 ) ;
}
