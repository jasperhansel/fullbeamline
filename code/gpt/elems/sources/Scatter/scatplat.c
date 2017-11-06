/* plate.c - Define a rectangular plate boundary */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

struct plate_info
{
  double a ;
  double b ;
} ;

static int plate_bound(gptpar *par, double t, double dt, gpttrajectory *trajectory, struct plate_info *info ) ;

void scatterplate_init(gptinit *init)
{
  struct plate_info *info ;
  int numarg ;

  gptbuildECS( init ) ;

  numarg = gptgetargnum(init) ;

  if( numarg!=2 )
    gpterror( "Syntax: %s(ECS,a,b)\n", gptgetname(init) ) ;

  info = (struct plate_info *)gptmalloc( sizeof(struct plate_info) ) ;

  info->a = gptgetargdouble(init,1)/2 ;
  info->b = gptgetargdouble(init,2)/2 ;

  gptaddboundaryelement( init, plate_bound, gptfree, 0, info ) ;
}


static int plate_bound(gptpar *par, double t, double dt, gpttrajectory *trajectory, struct plate_info *info )
{
  int i ;
  double a, b ;
  double lambda ;
  double P[3] ;

  a = info->a ;
  b = info->b ;

  /* See if the particle makes any chance of crossing the plate */
  if( trajectory->rstart[2] < 0.0 && trajectory->rend[2] < 0.0 ) return( 0 ) ;
  if( trajectory->rstart[2] > 0.0 && trajectory->rend[2] > 0.0 ) return( 0 ) ;
  if( trajectory->rstart[0] < -a && trajectory->rend[0] < -a ) return( 0 ) ;
  if( trajectory->rstart[0] >  a && trajectory->rend[0] >  a ) return( 0 ) ;
  if( trajectory->rstart[1] < -b && trajectory->rend[1] < -b ) return( 0 ) ;
  if( trajectory->rstart[1] >  b && trajectory->rend[1] >  b ) return( 0 ) ;

  /* Calculate intersection point (P=rstart+lambda*dr, see Mathematica document */
  lambda = -trajectory->rstart[2]/trajectory->dr[2] ;
  if( lambda > trajectory->lambda ) return( 0 ) ;
  for(i=0 ; i<3 ; i++) P[i] = trajectory->rstart[i] + lambda*trajectory->dr[i] ;

  /* Test if P is within plate dimensions */
  if( fabs(P[0]) >= a || fabs(P[1]) >= b ) return( 0 ) ;

  /* Calculate normal at intersection point, see Mathematica document */

  /* Store result */
  for(i=0 ; i<3 ; i++) trajectory->P[i] = P[i] ;
  trajectory->lambda = lambda ;
  trajectory->n[0] =  0.0 ;
  trajectory->n[1] =  0.0 ;
  trajectory->n[2] = -1.0 ;  

  return( 1 ) ;
}
