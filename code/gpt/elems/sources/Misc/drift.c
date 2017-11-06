/* drift.c - drift section */

/* BCM-General Particle Simulation: Drift section
 *
 * This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

struct drift_info
{
  double zlen ;
  double radius2 ;
} ;

static int drift_sim (gptpar *par,double t,struct drift_info *info) ;
static int drift_simr(gptpar *par,double t,struct drift_info *info) ;


void drift_init(gptinit *init)
{
  struct drift_info *info ;
  int numarg ;

  gptbuildECS( init ) ;

  numarg = gptgetargnum(init) ;

  if( numarg!=1 && numarg!=2 )
    gpterror( "Syntax: %s(ECS,L,[R])\n", gptgetname(init) ) ;

  info = (struct drift_info *)gptmalloc( sizeof(struct drift_info) ) ;

  info->zlen = gptgetargdouble(init,1)/2 ;

  if( numarg==2 )
  {
    info->radius2 = SQR(gptgetargdouble(init,2)) ;
    gptaddEBelement( init, drift_simr, gptfree, GPTELEM_LOCAL, info ) ;
  } else
    gptaddEBelement( init, drift_sim , gptfree, GPTELEM_LOCAL, info ) ;
}


static int drift_sim(gptpar *par,double t,struct drift_info *info)
{
  if( fabs(Z) > info->zlen ) return( 0 ) ;
  return( 1 ) ;
}


static int drift_simr(gptpar *par,double t,struct drift_info *info)
{
  if( fabs(Z) > info->zlen ) return( 0 ) ;
  if( SQR(X)+SQR(Y) <= info->radius2 ) return( 1 ) ;

  gptremoveparticle(par) ;
  return( 1 ) ;
}
