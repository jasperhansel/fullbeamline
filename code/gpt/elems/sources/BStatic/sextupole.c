/* sextupol.c - Sextupole lens */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

struct sextupole_info
{
  double G ;
  double zlen ;
} ;

static int sextupole_sim(gptpar *par,double t,struct sextupole_info *info) ;


void sextupole_init(gptinit *init)
{
  struct sextupole_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=2 )
    gpterror( "Syntax: %s(ECS,L,G)\n", gptgetname(init) ) ;

  info = (struct sextupole_info *)gptmalloc( sizeof(struct sextupole_info) ) ;

  info->zlen = gptgetargdouble(init,1)/2 ;
  info->G    = gptgetargdouble(init,2) ;

  gptaddEBelement( init, sextupole_sim, gptfree, GPTELEM_LOCAL, info ) ;
}


static int sextupole_sim(gptpar *par,double t,struct sextupole_info *info)
{
  if( fabs(Z) > info->zlen ) return( 0 ) ;

  BX = 2 * info->G * X * Y ;
  BY = info->G * (X*X-Y*Y) ;

  return( 1 ) ;
}
