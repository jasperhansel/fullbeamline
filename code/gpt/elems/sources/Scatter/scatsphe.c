/* sphere.c - Kill particle on a sphere and generate new reflected particles */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "elem.h"

#define EXTRAPOLATE 1e-6

struct sphere_info
{
  double R ;			/* Sphere radius */
  double a1, a2, adiff ;	/* Start, end and diffangle, 0 when not used */
} ;

static int sphere_bound(gptpar *par, double t, double dt, gpttrajectory *trajectory, struct sphere_info *info ) ;

void scattersphere_init(gptinit *init)
{
  struct sphere_info *info ;
  int numarg ;

  gptbuildECS( init ) ;

  numarg = gptgetargnum(init) ;

  if( numarg!=1 && numarg!=3 )
    gpterror( "Syntax: %s(ECS,R,[a1,a2])\n", gptgetname(init) ) ;

  info = (struct sphere_info *)gptmalloc( sizeof(struct sphere_info) ) ;

  info->R = gptgetargdouble(init,1) ;

  info->a1 = 0.0 ;
  info->a2 = 0.0 ;
  if( numarg==3 )
  {
    info->a1 = gptgetargdouble(init,2) ;
    info->a2 = gptgetargdouble(init,3) ;
    info->adiff = info->a2 - info->a1 ;
    if( info->adiff<-gpt_pi ) info->adiff += 2*gpt_pi ;
    if( info->adiff>+gpt_pi ) info->adiff -= 2*gpt_pi ;
    if( info->adiff<-0.999*gpt_pi || info->adiff>+0.999*gpt_pi ) info->adiff = gpt_pi ;
  }

  if( info->R<=0 )
    gpterror( "Sphere radius must be larger than zero\n" ) ;

  gptaddboundaryelement( init, sphere_bound, gptfree, 0, info ) ;
}


static int sphere_bound(gptpar *par, double t, double dt, gpttrajectory *trajectory, struct sphere_info *info )
{
  int i, count ;
  double R, Rstart2, Rend2 ;
  double *dr, result[2] ;
  double sx, sy, sz, a, b, c, lambda ;
  double P[3] ;
  double alpha ;

  R = info->R ;

  /* See if the particle makes any chance of crossing the sphere */
  if( trajectory->rstart[2] < -R && trajectory->rend[2] < -R ) return( 0 ) ;
  if( trajectory->rstart[2] >  R && trajectory->rend[2] >  R ) return( 0 ) ;
  if( trajectory->rstart[0] < -R && trajectory->rend[0] < -R ) return( 0 ) ;
  if( trajectory->rstart[0] >  R && trajectory->rend[0] >  R ) return( 0 ) ;
  if( trajectory->rstart[1] < -R && trajectory->rend[1] < -R ) return( 0 ) ;
  if( trajectory->rstart[1] >  R && trajectory->rend[1] >  R ) return( 0 ) ;

  Rstart2= gptVECSQR(trajectory->rstart) ;
  Rend2  = gptVECSQR(trajectory->rend) ;
  if( Rstart2 < R*R && Rend2 < R*R ) return( 0 ) ;

  /* Calculate trajectory direction dr=rend-rstart */
  dr = trajectory->dr ;

  /* Calculate intersection point (P=rstart+lambda*dr, see Mathematica document */
  sx = trajectory->rstart[0] ;
  sy = trajectory->rstart[1] ;
  sz = trajectory->rstart[2] ;
  a = (dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]) ;   /* lambda^2 */
  b = 2*( dr[0]*sx + dr[1]*sy + dr[2]*sz ) ;   /* lambda^1 */
  c = -R*R + sx*sx + sy*sy + sz*sz ;   /* lambda^0 */

  count = dblsolvequadratic( a, b, c, result ) ;
  if( count==0 ) return( 0 ) ;
  if( count==2 ) if( result[1]<0 || result[1]>1 ) count=1 ;
  if( count==2 ) if( result[0]<0 || result[0]>1 ) { result[0] = result[1] ; count=1 ; }
  if( count==2 ) if( result[0]>result[1] ) result[0]=result[1] ;

  lambda = result[0] ;
  if( lambda<0 || lambda>1.0+EXTRAPOLATE ) return( 0 ) ;
  if( lambda > trajectory->lambda ) return( 0 ) ;
  for(i=0 ; i<3 ; i++) P[i] = trajectory->rstart[i] + lambda*dr[i] ;

  /* Test angles */
  if( info->a1!=0.0 || info->a2!=0.0 )
  {
    alpha = atan2(sqrt(P[0]*P[0]+P[1]*P[1]),P[2]) - info->a1 ;
    if( alpha<-gpt_pi ) alpha += 2*gpt_pi ;
    if( alpha>+gpt_pi ) alpha -= 2*gpt_pi ;
    if( info->adiff * alpha < 0 ) return(0) ;
    if( fabs(alpha)>fabs(info->adiff) ) return(0) ;
  }

  /* Store result */
  for(i=0 ; i<3 ; i++) trajectory->P[i] = P[i] ;
  trajectory->lambda = lambda ;
  trajectory->n[0] =  P[0] ;
  trajectory->n[1] =  P[1] ;
  trajectory->n[2] =  P[2] ;

  return( 1 ) ;
}
