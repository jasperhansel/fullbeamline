/* Pulsar Physics: Sector magnet with non-normal pole faces and fringe fields */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 *
 * With special thanks to Bruno Muratory from Daresbury Laboratory
 * for providing us with the fringe field expressions
 */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "elem.h"

#define NGAP 10 /* Transient region is NGAP/b1 */
#define FUZZFACTOR (1024*DBL_EPSILON)

struct sectormagnet_info
{
  double n[3] ;  /* Vector normal to center plane */
  double m[3] ;  /* Vector normal to pole-face */
  double B ;     /* Magnetic field strength */
  double f[3] ;  /* Magnetic field direction, normal to plane of fromaxis and toaxis */
  double l ;	 /* Distance from axis intersection to edge */
  double dl ;    /* Offset for edge */
  double b1, b2 ;/* Parameters of Enge fringe field */
  struct axis *otheraxis ;
} ;

static int sectormagnet_sim(gptpar *par,double t,struct sectormagnet_info *info) ;
static int sectormagnet_simple_sim(gptpar *par,double t,struct sectormagnet_info *info) ;

void sectormagnet_init(gptinit *init)
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
  double m[3] ;       /* Vector normal to pole-face */
  double B ;	      /* Specified B field */
  double phi_in, phi_out ;	/* in- and output angles in bend plane */
  double dl ;         /* Offset for edge (difference between real and effective length) */
  double b1, b2 ;     /* Parameters of Enge fringe field */ 

  char *fromaxisname, *toaxisname ;
  struct axis *fromaxis, *toaxis ;
  double de, dpq, epq, Pd, Pe ;
  double divf ;
  int argc ;
  unsigned int i, j ;
  struct sectormagnet_info *info ;

  argc = gptgetargnum(init) ;
  if( argc!=4 && argc!=9)
    gpterror( "Syntax: %s(fromaxis,toaxis,R,Bfield,[phi_in,phi_out,dl,b1,b2])\n", gptgetname(init) ) ;

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

/* Get angles of pole-faces */
  phi_in = phi_out = 0.0 ;
  dl = 0 ;
  b1 = b2 = 0.0 ;
  if( argc==9 )
  {
    phi_in  = gptgetargdouble(init,5) ;
    phi_out = gptgetargdouble(init,6) ;
    dl      = gptgetargdouble(init,7) ;
    b1      = gptgetargdouble(init,8) ;
    b2      = gptgetargdouble(init,9) ;
  }

  if( b1<0 )
    gpterror( "%s: b1 must be nonnegative\n", gptgetname(init) ) ;

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
  info = (struct sectormagnet_info *)gptmalloc( sizeof(struct sectormagnet_info) ) ;

  init->e.o[2] = t ;
  init->paxis = fromaxis ;
  info->otheraxis = toaxis ;
  m[0] = cos(phi_in)*d[0] - sin(phi_in)*(d[1]*f[2]-d[2]*f[1]) ; /* + tan(theta_in)*f[0] ; */
  m[1] = cos(phi_in)*d[1] - sin(phi_in)*(d[2]*f[0]-d[0]*f[2]) ; /* + tan(theta_in)*f[1] ; */
  m[2] = cos(phi_in)*d[2] - sin(phi_in)*(d[0]*f[1]-d[1]*f[0]) ; /* + tan(theta_in)*f[2] ; */
  BtoUCS( &fromaxis->a, n, info->n ) ;
  BtoUCS( &fromaxis->a, f, info->f ) ;
  BtoUCS( &fromaxis->a, m, info->m ) ;
  info->B  = B ;
  info->l  = l ;
  info->dl = dl ;
  info->b1 = b1 ;
  info->b2 = b2 ;

  if( info->b1==0 )
    gptaddEBelement( init, sectormagnet_simple_sim, gptfree, GPTELEM_LOCAL, info ) ;
  else
    gptaddEBelement( init, sectormagnet_sim, gptfree, GPTELEM_GLOBAL, info ) ;


/* Add to toaxis */
  info = (struct sectormagnet_info *)gptmalloc( sizeof(struct sectormagnet_info) ) ;

  init->e.o[2] = u ;
  init->paxis = toaxis ;
  info->otheraxis = fromaxis ;
  m[0] = cos(phi_out)*e[0] + sin(phi_out)*(e[1]*f[2]-e[2]*f[1]) ; /* - tan(theta_in)*f[0] ; */
  m[1] = cos(phi_out)*e[1] + sin(phi_out)*(e[2]*f[0]-e[0]*f[2]) ; /* - tan(theta_in)*f[1] ; */
  m[2] = cos(phi_out)*e[2] + sin(phi_out)*(e[0]*f[1]-e[1]*f[0]) ; /* - tan(theta_in)*f[2] ; */
  BtoUCS( &toaxis->a, n, info->n ) ; for(i=0 ; i<3 ; i++) info->n[i] *= -1 ;
  BtoUCS( &toaxis->a, f, info->f ) ;
  BtoUCS( &toaxis->a, m, info->m ) ; for(i=0 ; i<3 ; i++) info->m[i] *= -1 ;
  info->B  = B ;
  info->l  = -l ;
  info->dl = dl ;
  info->b1 = b1 ;
  info->b2 = b2 ;

  if( info->b1==0 )
    gptaddEBelement( init, sectormagnet_simple_sim, gptfree, GPTELEM_LOCAL, info ) ;
  else
    gptaddEBelement( init, sectormagnet_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}


static int sectormagnet_sim(gptpar *par,double t,struct sectormagnet_info *info)
{
  double d ;    /* Distance from edge (-Z)*/
  double e ;    /* Distance in magnetic field direction (Y, when particle bends in XZ) */
  double Bf ;	/* Magnetic field in e direction */
  double Bm ;	/* Magnetic field in d direction */
  double P[3] ; /* Intersection point between edge and axis */
  double f, h ;
  double Bn ;

  P[0] = X ;
  P[1] = Y ;
  P[2] = Z+info->l ;

  d = -gptVECINP(P,info->m) - info->dl;
  if( d*info->b1 > NGAP ) return( 0 ) ;

  e = gptVECINP(par->r,info->f) ;
  if( fabs(e)*info->b1>=gpt_pi && d<=0 ) { gptremoveparticle(par) ; return(1) ; }

  f = info->b1*d + info->b2*(d*d-e*e) ;
  h = e*(info->b1+2*info->b2*d) ;

  Bn = info->B/(1+2*exp(f)*cos(h)+exp(2*f)) ;

  Bf = Bn*(1+exp(f)*cos(h)) ;
  Bm = Bn*exp(f)*sin(h) ;

  BX = Bf*info->f[0] + Bm*info->m[0] ;
  BY = Bf*info->f[1] + Bm*info->m[1] ;
  BZ = Bf*info->f[2] + Bm*info->m[2] ;

  if( gptVECINP(par->r, info->n) < 0 ) par->newaxis = info->otheraxis ;

  return( 1 ) ;
}


static int sectormagnet_simple_sim(gptpar *par,double t,struct sectormagnet_info *info)
{
  double d ;    /* Distance from edge */
  double P[3] ; /* Intersection point between edge and axis */

  P[0] = X ;
  P[1] = Y ;
  P[2] = Z+info->l ;

  d = -gptVECINP(P,info->m) ;
  if( d > 0.0 ) return( 0 ) ;

  BX = info->B*info->f[0] ;
  BY = info->B*info->f[1] ;
  BZ = info->B*info->f[2] ;

  if( gptVECINP(par->r, info->n) < 0 ) par->newaxis = info->otheraxis ;

  return( 1 ) ;
}
