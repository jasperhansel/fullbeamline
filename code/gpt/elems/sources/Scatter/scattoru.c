/* torus.c - Kill particle on a torus and generate new reflected particles */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "elem.h"

#define EXTRAPOLATE 1e-6

struct torus_info
{
  double R1, R2 ;		/* Outer and inner radius resp. */
  double a1, a2, adiff ;	/* Start, end and diffangle, 0 when not used */
} ;

static int torus_bound(gptpar *par, double t, double dt, gpttrajectory *trajectory, struct torus_info *info ) ;
static int dblsolvecube4simple(double a, double b, double c, double d, double e, double *result ) ;

void scattertorus_init(gptinit *init)
{
  struct torus_info *info ;
  int numarg ;

  gptbuildECS( init ) ;

  numarg = gptgetargnum(init) ;

  if( numarg!=2 && numarg!=4 )
    gpterror( "Syntax: %s(ECS,Rout,Rin,a1,a2)\n", gptgetname(init) ) ;

  info = (struct torus_info *)gptmalloc( sizeof(struct torus_info) ) ;

  info->R1 = gptgetargdouble(init,1) ;
  info->R2 = gptgetargdouble(init,2) ;

  info->a1 = 0.0 ;
  info->a2 = 0.0 ;
  if( numarg==4 )
  {
    info->a1 = gptgetargdouble(init,3) ;
    info->a2 = gptgetargdouble(init,4) ;
    info->adiff = info->a2 - info->a1 ;
    if( info->adiff<-gpt_pi ) info->adiff += 2*gpt_pi ;
    if( info->adiff>+gpt_pi ) info->adiff -= 2*gpt_pi ;
    if( info->adiff<-0.999*gpt_pi || info->adiff>+0.999*gpt_pi ) info->adiff = gpt_pi ;
  }

  if( info->R1<=0 || info->R2<=0 )
    gpterror( "Torus dimensions must be larger than zero\n" ) ;
#if 0
  if( info->R1<info->R2 )
    gpterror( "Outer radius of torus must be larger than inner radius\n" ) ;
  if( info->R1==info->R2 )
    gpterror( "Outer radius of torus cannot be equal to inner radius\n" ) ;
#endif

  gptaddboundaryelement( init, torus_bound, gptfree, 0, info ) ;
}


static int torus_bound(gptpar *par, double t, double dt, gpttrajectory *trajectory, struct torus_info *info )
{
  int i, count ;
  double R1, R2, R1plR2, Rstart2, Rend2 ;
  double *dr, result[4] ;
  double sx, sy, sz, a, b, c, d, e, lambda ;
  double P[3] ;
  double alpha ;

  R1 = info->R1 ;
  R2 = info->R2 ;
  R1plR2 = R1+R2 ;

  /* See if the particle makes any chance of crossing the torus */
  if( trajectory->rstart[2] < -R2 && trajectory->rend[2] < -R2 ) return( 0 ) ;
  if( trajectory->rstart[2] >  R2 && trajectory->rend[2] >  R2 ) return( 0 ) ;
  if( trajectory->rstart[0] < -R1plR2 && trajectory->rend[0] < -R1plR2 ) return( 0 ) ;
  if( trajectory->rstart[0] >  R1plR2 && trajectory->rend[0] >  R1plR2 ) return( 0 ) ;
  if( trajectory->rstart[1] < -R1plR2 && trajectory->rend[1] < -R1plR2 ) return( 0 ) ;
  if( trajectory->rstart[1] >  R1plR2 && trajectory->rend[1] >  R1plR2 ) return( 0 ) ;

  Rstart2= trajectory->rstart[0]*trajectory->rstart[0]+trajectory->rstart[1]*trajectory->rstart[1] ;
  Rend2  = trajectory->rend[0]*trajectory->rend[0]+trajectory->rend[1]*trajectory->rend[1] ;
  if( Rstart2 < (R1-R2)*(R1-R2) && Rend2 < (R1-R2)*(R1-R2) ) return( 0 ) ;

  /* Calculate trajectory direction dr=rend-rstart */
  dr = trajectory->dr ;

  /* Calculate intersection point (P=rstart+lambda*dr, see Mathematica document */
  sx = trajectory->rstart[0] ;
  sy = trajectory->rstart[1] ;
  sz = trajectory->rstart[2] ;
  a = (dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2])*(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]) ;   /* lambda^4 */
  b = 4*(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2])*( dr[0]*sx + dr[1]*sy + dr[2]*sz ) ;   /* lambda^3 */
  c = 2*(-R1*R1*(dr[0]*dr[0] + dr[1]*dr[1] - dr[2]*dr[2])-R2*R2*(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]) +
      4*(dr[0]*dr[1]*sx*sy + dr[0]*dr[2]*sx*sz + dr[1]*dr[2]*sy*sz) + 
      (3*dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2])*sx*sx + 
      (dr[0]*dr[0] + 3*dr[1]*dr[1] + dr[2]*dr[2])*sy*sy +
      (dr[0]*dr[0] + dr[1]*dr[1] + 3*dr[2]*dr[2])*sz*sz ) ;   /* lambda^2 */
  d = 4*(-R1*R1*(dr[0]*sx+dr[1]*sy-dr[2]*sz) + (dr[0]*sx+dr[1]*sy+dr[2]*sz)*(sx*sx+sy*sy+sz*sz-R2*R2)) ;   /* lambda^1 */
  e = R1*R1*R1*R1-2*R1*R1*(sx*sx+sy*sy-sz*sz+R2*R2) + (sx*sx+sy*sy+sz*sz-R2*R2)*(sx*sx+sy*sy+sz*sz-R2*R2) ;   /* lambda^0 */

  count = dblsolvecube4simple( a, b, c, d, e, result ) ;
  if( count!=1 ) return( 0 ) ;

  lambda = result[0] ;
  if( lambda<0 || lambda>1.0+EXTRAPOLATE ) return( 0 ) ;
  if( lambda > trajectory->lambda ) return( 0 ) ;
  for(i=0 ; i<3 ; i++) P[i] = trajectory->rstart[i] + lambda*dr[i] ;

  /* Test angles */
  if( info->a1!=0.0 || info->a2!=0.0 )
  {
    alpha = atan2(sqrt(P[0]*P[0]+P[1]*P[1])-R1,P[2]) - info->a1 ;
    if( alpha<-gpt_pi ) alpha += 2*gpt_pi ;
    if( alpha>+gpt_pi ) alpha -= 2*gpt_pi ;
    if( info->adiff * alpha < 0 ) return(0) ;
    if( fabs(alpha)>fabs(info->adiff) ) return(0) ;
  }

  /* Store result */
  for(i=0 ; i<3 ; i++) trajectory->P[i] = P[i] ;
  trajectory->lambda = lambda ;
  trajectory->n[0] =  P[0]*(-R1*R1-R2*R2+P[0]*P[0]+P[1]*P[1]+P[2]*P[2]) ;
  trajectory->n[1] =  P[1]*(-R1*R1-R2*R2+P[0]*P[0]+P[1]*P[1]+P[2]*P[2]) ;
  trajectory->n[2] =  P[2]*( R1*R1-R2*R2+P[0]*P[0]+P[1]*P[1]+P[2]*P[2]) ;  

  return( 1 ) ;
}


/* Solve ax^4+bx^3+cx^2+dx+e==0 for the smallest real solution bewteen 0 and 1
 *
 * Modified from Muller's Method in NR
 * Bugs: Could accidentally find a larger than the smallest root
 */

#define PREC4 1e-8
#define MAXIT4 16

static double realpoly4(double x, double a, double b, double c, double d, double e)
{
  return( a*x*x*x*x + b*x*x*x + c*x*x + d*x + e) ;
}

static int dblsolvecube4simple(double a, double b, double c, double d, double e, double *result )
{
  double x[3], P[3] ;
  int i ;

  x[0] = 0.0 ; P[0] = realpoly4(x[0],a,b,c,d,e) ;
  x[1] = 0.5 ; P[1] = realpoly4(x[1],a,b,c,d,e) ;
  x[2] = 1.0 ; P[2] = realpoly4(x[2],a,b,c,d,e) ;

  for(i=0 ; i<MAXIT4 ; i++)
  {
    double q, A, B, C ;
    double det ;

    if( fabs(x[0]-x[1])<PREC4 )
    {
      *result = x[0] ;
      if( realpoly4(x[0]+PREC4,a,b,c,d,e)*realpoly4(x[0]-PREC4,a,b,c,d,e) < 0 ) return(1) ;
      return(0) ;
    }
    if( x[0] < -1 || x[0] > 2 ) return(0) ;

    q = (x[0]-x[1])/(x[1]-x[2]) ;
    A = q*P[0] - q*(1+q)*P[1] + q*q*P[2] ;
    B = (2*q+1)*P[0] - (1+q)*(1+q)*P[1] + q*q*P[2] ;
    C = (1+q)*P[0] ;

    x[2] = x[1] ; x[1] = x[0] ;
    P[2] = P[1] ; P[1] = P[0] ;

    det = B*B-4*A*C ;
    if( B==0 && det==0 ) /* Stop and test the root */
      x[0] = x[1] ;
    else  if( det<0 )  /* Stay on real axis */
      x[0] = x[1] - (x[1]-x[2])*2*C*B/(B*B-det) ;
    else
      x[0] = x[1] - (x[1]-x[2])*2*C/(B>0 ? B+sqrt(det) : B-sqrt(det)) ;

    P[0] = realpoly4(x[0],a,b,c,d,e) ;
  }

  return( 1 ) ;
}
