/* cone.c - Kill particle on a cone and generate new reflected particles */

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

struct cone_info
{
  double r1, z1 ;
  double r2, z2 ;
  double rmin2, rmax ;
  double zmin, zmax ;
} ;

static int cone_bound(gptpar *par, double t, double dt, gpttrajectory *trajectory, struct cone_info *info ) ;

void scattercone_init(gptinit *init)
{
  struct cone_info *info ;
  int numarg ;

  gptbuildECS( init ) ;

  numarg = gptgetargnum(init) ;

  if( numarg!=4 )
    gpterror( "Syntax: %s(ECS,z1,R1,z2,R2)\n", gptgetname(init) ) ;

  info = (struct cone_info *)gptmalloc( sizeof(struct cone_info) ) ;

  info->z1 = gptgetargdouble(init,1) ;
  info->r1 = gptgetargdouble(init,2) ;
  info->z2 = gptgetargdouble(init,3) ;
  info->r2 = gptgetargdouble(init,4) ;

  if( info->z1==info->z2 )
    gpterror( "z1 and z2 equal to %lg. Use iris element.\n", info->z1 ) ;

  if( info->r1<0 || info->r2<0 )
    gpterror( "Cone dimensions must be nonnegative\n" ) ;

  if( info->r1==info->r2 )
    gpterror( "r1 and r2 equal to %lg. Use pipe element.\n", info->r1 ) ;

  info->rmin2 = dblMIN(info->r1,info->r2)*dblMIN(info->r1,info->r2) ;
  info->rmax  = dblMAX(info->r1,info->r2) ;
  info->zmin = dblMIN(info->z1,info->z2) ;
  info->zmax = dblMAX(info->z1,info->z2) ;

  gptaddboundaryelement( init, cone_bound, gptfree, 0, info ) ;
}


static int cone_bound(gptpar *par, double t, double dt, gpttrajectory *trajectory, struct cone_info *info )
{
  int i, count ;
  double r1, z1, r2, z2, rmax, Rstart2, Rend2, PR2 ;
  double *dr, result[2] ;
  double sx, sy, sz, alpha, a, b, c, lambda ;
  double P[3] ;

  r1 = info->r1 ;
  z1 = info->z1 ;
  r2 = info->r2 ;
  z2 = info->z2 ;
  rmax = info->rmax ;

  /* See if the particle makes any chance of crossing the cone */
  if( trajectory->rstart[2] < info->zmin && trajectory->rend[2] < info->zmin ) return( 0 ) ;
  if( trajectory->rstart[2] > info->zmax && trajectory->rend[2] > info->zmax ) return( 0 ) ;
  if( trajectory->rstart[0] < -rmax && trajectory->rend[0] < -rmax ) return( 0 ) ;
  if( trajectory->rstart[0] >  rmax && trajectory->rend[0] >  rmax ) return( 0 ) ;
  if( trajectory->rstart[1] < -rmax && trajectory->rend[1] < -rmax ) return( 0 ) ;
  if( trajectory->rstart[1] >  rmax && trajectory->rend[1] >  rmax ) return( 0 ) ;

  Rstart2= trajectory->rstart[0]*trajectory->rstart[0]+trajectory->rstart[1]*trajectory->rstart[1] ;
  Rend2  = trajectory->rend[0]*trajectory->rend[0]+trajectory->rend[1]*trajectory->rend[1] ;
  if( Rstart2 < info->rmin2 && Rend2 < info->rmin2 ) return( 0 ) ;

  /* Calculate trajectory direction dr=rend-rstart */
  dr = trajectory->dr ;

  /* Calculate intersection point (P=rstart+lambda*dr, see Mathematica document */
  sx = trajectory->rstart[0] ;
  sy = trajectory->rstart[1] ;
  sz = trajectory->rstart[2] ;
  alpha  = (r2-r1)/(z2-z1) ;
  a = dr[0]*dr[0] + dr[1]*dr[1] - alpha*alpha*dr[2]*dr[2] ;
  b = 2*( dr[0]*sx + dr[1]*sy - alpha*dr[2]*(alpha*sz-alpha*z1+r1) ) ;
  c = sx*sx + sy*sy - (alpha*sz-alpha*z1+r1)*(alpha*sz-alpha*z1+r1) ;

  count = dblsolvequadratic( a, b, c, result ) ;
  if( count==0 ) return( 0 ) ;
  if( count==2 ) if( result[1]<0 || result[1]>1 ) count=1 ;
  if( count==2 ) if( result[0]<0 || result[0]>1 ) { result[0] = result[1] ; count=1 ; }
  if( count==2 ) if( result[0]>result[1] ) result[0]=result[1] ;

  lambda = result[0] ;
  if( lambda<0 || lambda>1.0+EXTRAPOLATE ) return( 0 ) ;
  if( lambda > trajectory->lambda ) return( 0 ) ;
  for(i=0 ; i<3 ; i++) P[i] = trajectory->rstart[i] + lambda*dr[i] ;

  /* Test if P is within cone dimensions */
  if( P[2] < info->zmin || P[2] > info->zmax ) return( 0 ) ;
  PR2 = P[0]*P[0] + P[1]*P[1] ;
  if( PR2 < info->rmin2 || PR2 > rmax*rmax ) return( 0 ) ;

  /* Store results */
  for(i=0 ; i<3 ; i++) trajectory->P[i] = P[i] ;
  trajectory->lambda = lambda ;
  trajectory->n[0] =  P[0] ;
  trajectory->n[1] =  P[1] ;
  trajectory->n[2] =  alpha*(r2*(P[2]-z1) + r1*(z2-P[2]))/(z1-z2) ;  

  return( 1 ) ;
}
