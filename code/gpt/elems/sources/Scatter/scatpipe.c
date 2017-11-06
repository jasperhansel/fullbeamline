/* pipe.c - Kill particle on a pipe and generate new reflected particles */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "elem.h"

#define dblMIN(a,b) ((a)>(b)?(b):(a))
#define dblMAX(a,b) ((a)>(b)?(a):(b))
#define EXTRAPOLATE 1e-6

struct pipe_info
{
  double R ;
  double zmin, zmax ;
} ;

static int pipe_bound(gptpar *par, double t, double dt, gpttrajectory *trajectory, struct pipe_info *info ) ;

void scatterpipe_init(gptinit *init)
{
  double z1, z2 ;
  struct pipe_info *info ;
  int numarg ;

  gptbuildECS( init ) ;

  numarg = gptgetargnum(init) ;

  if( numarg!=3 )
    gpterror( "Syntax: %s(ECS,z1,z2,R)\n", gptgetname(init) ) ;

  info = (struct pipe_info *)gptmalloc( sizeof(struct pipe_info) ) ;

  z1 = gptgetargdouble(init,1) ;
  z2 = gptgetargdouble(init,2) ;
  info->R  = gptgetargdouble(init,3) ;

  if( info->R==0 )
    return ;
  if( z1==z2 )
    gpterror( "z1 and z2 equal to %lg\n", z1 ) ;
  if( info->R<0 )
    gpterror( "Pipe dimensions must be larger than zero\n" ) ;

  info->zmin = dblMIN(z1,z2) ;
  info->zmax = dblMAX(z1,z2) ;

  gptaddboundaryelement( init, pipe_bound, gptfree, 0, info ) ;
}


static int pipe_bound(gptpar *par, double t, double dt, gpttrajectory *trajectory, struct pipe_info *info )
{
  int i, count ;
  double R, Rstart2, Rend2 ;
  double *dr, result[2] ;
  double sx, sy, sz, a, b, c, lambda ;
  double P[3] ;

  R  = info->R ;

  /* See if the particle makes any chance of crossing the pipe */
  if( trajectory->rstart[2] < info->zmin && trajectory->rend[2] < info->zmin ) return( 0 ) ;
  if( trajectory->rstart[2] > info->zmax && trajectory->rend[2] > info->zmax ) return( 0 ) ;
  if( trajectory->rstart[0] < -R && trajectory->rend[0] < -R ) return( 0 ) ;
  if( trajectory->rstart[0] >  R && trajectory->rend[0] >  R ) return( 0 ) ;
  if( trajectory->rstart[1] < -R && trajectory->rend[1] < -R ) return( 0 ) ;
  if( trajectory->rstart[1] >  R && trajectory->rend[1] >  R ) return( 0 ) ;

  Rstart2= gptSQR(trajectory->rstart[0])+gptSQR(trajectory->rstart[1]) ;
  Rend2  = gptSQR(trajectory->rend[0])  +gptSQR(trajectory->rend[1])   ;
  if( Rstart2 < R*R && Rend2 < R*R ) return( 0 ) ; 

  /* Calculate trajectory direction dr=rend-rstart */
  dr = trajectory->dr ;

  /* Calculate intersection point (P=rstart+lambda*dr, see Mathematica document */
  sx = trajectory->rstart[0] ;
  sy = trajectory->rstart[1] ;
  sz = trajectory->rstart[2] ;
  a = dr[0]*dr[0] + dr[1]*dr[1] ;
  b = 2*( dr[0]*sx + dr[1]*sy ) ;
  c = sx*sx + sy*sy - R*R ;

  count = dblsolvequadratic( a, b, c, result ) ;
  if( count==0 ) return( 0 ) ;
  if( count==2 ) if( result[1]<0 || result[1]>1 ) count=1 ;
  if( count==2 ) if( result[0]<0 || result[0]>1 ) { result[0] = result[1] ; count=1 ; }
  if( count==2 ) if( result[0]>result[1] ) result[0]=result[1] ;

  lambda = result[0] ;
  if( lambda<0 || lambda>1.0+EXTRAPOLATE ) return( 0 ) ;
  if( lambda > trajectory->lambda ) return( 0 ) ;
  for(i=0 ; i<3 ; i++) P[i] = trajectory->rstart[i] + lambda*dr[i] ;

  /* Test if P is within pipe dimensions */
  if( P[2] < info->zmin || P[2] > info->zmax ) return( 0 ) ;

  /* Store result */
  for(i=0 ; i<3 ; i++) trajectory->P[i] = P[i] ;
  trajectory->lambda = lambda ;
  trajectory->n[0] =  P[0] ;
  trajectory->n[1] =  P[1] ;
  trajectory->n[2] =  0 ;  

  return( 1 ) ;
}
