/* bend.c - bending magnet */

/* BCM-General Particle Simulation: Bending magnet
 *
 * This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "elem.h"

#define FUZZFACTOR (1024*DBL_EPSILON)

struct bend_info
{
  double n[3] ;
  double fB[3] ;
  double l ;
  struct axis *otheraxis ;
} ;

static int bend_sim1(gptpar *par,double t,struct bend_info *info) ;
static int bend_sim2(gptpar *par,double t,struct bend_info *info) ;

void bend_init(gptinit *init)
{
  double P[3], Q[3] ; /* Nearest poits on axis, P==Q for intersection */
  double *po, *qo ;   /* Origins */
  double d[3], e[3] ; /* Z-directions */
  double t, u ;       /* P=p+td, Q=q+ue */
  double a, b ;       /* C=ad+bd */
  double f[3] ;       /* f=d x e, LEN(f)=1 */
  double l ;	      /* Half of the length */
  double r ;	      /* Radius specified by user */
  double PC[3] ;      /* P - Center of circle */
  double n[3] ;       /* Vector normal to center-of-magnet plane */
  double B ;	      /* Specified B field */

  char *fromaxisname, *toaxisname ;
  struct axis *fromaxis, *toaxis ;
  double de, dpq, epq, Pd, Pe ;
  double divf ;
  unsigned int i, j ;
  struct bend_info *info ;

  if( gptgetargnum(init)!=4 )
    gpterror( "Syntax: %s(fromaxis, toaxis, R, Bfield)\n", gptgetname(init) ) ;

/* Get axisnames */
  fromaxisname=gptgetargstring(init,1) ;
  toaxisname  =gptgetargstring(init,2) ;
  if( (fromaxis=getaxis(fromaxisname))==NULL )
    gpterror( "Specified axis \"%s\" not found\n", fromaxisname ) ;
  if( (  toaxis=getaxis(  toaxisname))==NULL )
    gpterror( "Specified axis \"%s\" not found\n", toaxisname ) ;

/* Fill coefficients for local use. */
  r = gptgetargdouble(init,3) ;
  B = gptgetargdouble(init,4) ;
  po = fromaxis->a.o ;
  qo =	 toaxis->a.o ;
  for(i=0 ; i<3 ; i++) d[i] = fromaxis->a.m[i][2] ;
  for(i=0 ; i<3 ; i++) e[i] =	toaxis->a.m[i][2] ;

/* Calculate P and Q */
  de = INP(d,e) ;
  if( (divf=1.0-SQR(de)) < FUZZFACTOR )
    gpterror( "Axes point in the same direction\n" ) ;
  dpq = INP(d,po)-INP(d,qo) ;
  epq = INP(e,po)-INP(e,qo) ;
  t = (de*epq - dpq)/divf ;
  u = (epq - de*dpq)/divf ;
  for(i=0 ; i<3 ; i++) P[i] = po[i]+t*d[i] ;
  for(i=0 ; i<3 ; i++) Q[i] = qo[i]+u*e[i] ;

/* Check intersection */
  if( sqrt( SQR(P[0]-Q[0]) + SQR(P[1]-Q[1]) + SQR(P[2]-Q[2]) ) > FUZZFACTOR )
    gpterror( "Axes do not intersect. Nearest coordinates in WCS:\n"
    "%e %e %e on axis %s\n%e %e %e on axis %s\n",
    P[0], P[1], P[2], fromaxisname, Q[0], Q[1], Q[2], toaxisname ) ;

/* Calculate f and l */
  f[0] = d[1]*e[2] - d[2]*e[1] ;
  f[1] = d[2]*e[0] - d[0]*e[2] ;
  f[2] = d[0]*e[1] - d[1]*e[0] ;
  divf = sqrt( INP(f,f) ) ;
  for(i=0 ; i<3 ; i++) f[i]/=divf ;
  l = r*divf/(1.0+de) ;

/* Calculate PC */
  Pd = INP(P,d) ;
  Pe = INP(P,e) ;
  a = (Pd-l-de*(Pe+l))/(1.0-de*de) ;
  b = (Pe+l-de*(Pd-l))/(1.0-de*de) ;
  for(i=0 ; i<3 ; i++) PC[i] = P[i] - a*d[i] - b*e[i] ;

/* calculate n */
  n[0] = PC[1]*f[2] - PC[2]*f[1] ;
  n[1] = PC[2]*f[0] - PC[0]*f[2] ;
  n[2] = PC[0]*f[1] - PC[1]*f[0] ;

/* Fill in information for both bends */
  init->e.o[0] = init->e.o[1] = 0.0 ;
  for(i=0 ; i<3 ; i++) for(j=0 ; j<3 ; j++)
    if( i==j ) init->e.m[i][j] = 1.0 ; else init->e.m[i][j] = 0.0 ;
  init->e.mid = 1 ;

/* Add to fromaxis */
  info = (struct bend_info *)gptmalloc( sizeof(struct bend_info) ) ;

  init->e.o[2] = t ;
  init->paxis = fromaxis ;
  info->otheraxis = toaxis ;
  BtoUCS( &fromaxis->a, n, info->n ) ;
  BtoUCS( &fromaxis->a, f, info->fB ) ;
  for(i=0 ; i<3 ; i++) info->fB[i] *= B ;
  info->l = l ;

  gptaddEBelement( init, bend_sim1, gptfree, GPTELEM_LOCAL, info ) ;

/* Add to toaxis */
  info = (struct bend_info *)gptmalloc( sizeof(struct bend_info) ) ;

  init->e.o[2] = u ;
  init->paxis = toaxis ;
  info->otheraxis = fromaxis ;
  BtoUCS( &toaxis->a, n, info->n ) ;
  BtoUCS( &toaxis->a, f, info->fB ) ;
  for(i=0 ; i<3 ; i++) info->fB[i] *= B ;
  info->l = l ;

  gptaddEBelement( init, bend_sim2, gptfree, GPTELEM_LOCAL, info ) ;
}

static int bend_sim1(gptpar *par,double t,struct bend_info *info)
{
  if( Z < -info->l ) return( 0 ) ;

  BX = info->fB[0] ;
  BY = info->fB[1] ;
  BZ = info->fB[2] ;

  if( INP(par->r, info->n) < 0 ) par->newaxis = info->otheraxis ;

  return( 1 ) ;
}


static int bend_sim2(gptpar *par,double t,struct bend_info *info)
{
  if( Z > info->l ) return( 0 ) ;

  BX = info->fB[0] ;
  BY = info->fB[1] ;
  BZ = info->fB[2] ;

  if( INP(par->r, info->n) > 0 ) par->newaxis = info->otheraxis ;

  return( 1 ) ;
}

