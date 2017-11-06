/* bz.c - tabulated Bz element using splines */

/* BCM-General Particle Simulation: Tabulated Bz-file, using spline
 *
 * This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

struct bz_info
{
  unsigned int n ;
  double *z ;
  double *Bz ;
  double *Bz2 ;
} ;

static int bz_sim(gptpar *par,double t,struct bz_info *info) ;
static void bz_exit(struct bz_info *info) ;

void bz_init(gptinit *init)
{
  struct bz_info *info ;
  double current ;
  unsigned int i ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=2 )
    gpterror( "Syntax: %s(ECS, filename, I)\n", gptgetname(init) ) ;

  info = (struct bz_info *)gptmalloc( sizeof(struct bz_info) ) ;

  current=gptgetargdouble(init,2) ;

  readtabfile( init, gptgetargstring(init,1), &info->z, &info->Bz, &info->n ) ;

  info->Bz2 = (double *)gptmalloc(info->n*sizeof(double)) ;

  for(i=0 ; i<info->n ; i++) info->Bz[i] *= current ;

  if( spline( info->z, info->Bz, info->Bz2, info->n, NATSPLINE, NATSPLINE ) )
    gpterror( "Error calculating spline\n" ) ;

  gptaddEBelement( init, bz_sim, bz_exit, GPTELEM_GLOBAL, info ) ;
}


static int bz_sim(gptpar *par,double t,struct bz_info *info)
{
  double y, dy ;

  if( splintd( info->z, info->Bz, info->Bz2, info->n, Z,
    &y, &dy ) ) return( 0 ) ;

  BX = -dy*X * 0.5 ;
  BY = -dy*Y * 0.5 ;
  BZ = y ;

  return( 1 ) ;
}


static void bz_exit( struct bz_info *info )
{
  gptfree( info->z ) ;
  gptfree( info->Bz ) ;
  gptfree( info->Bz2 ) ;
  gptfree( info ) ;
}
