/* linecurent.c - Line with homogeneous current */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "elem.h"

#define ffac (gpt_mu0/(4*gpt_pi))

struct linecurrent_info
{
  double mfac ;
  double L ;
} ;

static int buildxyfromz(gpttransform *t) ;
static int linecurrent_sim(gptpar *par,double t,struct linecurrent_info *info) ;


void linecurrent_init(gptinit *init)
{
  struct linecurrent_info *info ;
  double r1[3], r2[3], r2mr1[3] ;
  int i ;
  gpttransform N, MN ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=7 )
    gpterror( "Syntax: %s(ECS,x1,y1,z1,x2,y2,z2,I)\n", gptgetname(init) ) ;

  info = (struct linecurrent_info *)gptmalloc( sizeof(struct linecurrent_info) ) ;

  r1[0]      = gptgetargdouble(init,1) ;
  r1[1]      = gptgetargdouble(init,2) ;
  r1[2]      = gptgetargdouble(init,3) ;
  r2[0]      = gptgetargdouble(init,4) ;
  r2[1]      = gptgetargdouble(init,5) ;
  r2[2]      = gptgetargdouble(init,6) ;
  info->mfac = ffac*gptgetargdouble(init,7) ;

  /* Optimalisation: Return when no current */
  if( info->mfac==0.0 ) return ;

  /* Calculate length of line */
  for( i=0 ; i<3 ; i++ ) r2mr1[i] = r2[i]-r1[i] ;
  info->L = gptVECLEN(r2mr1)/2 ;

  /* Calculate coordinate transformation */
  for( i=0 ; i<3 ; i++ ) N.o[i] = (r1[i]+r2[i])/2 ;
  for( i=0 ; i<3 ; i++ ) N.m[i][2] = r2mr1[i] ;
  if( buildxyfromz(&N)!=0)
    gpterror( "%s: Start and end points must be different\n", gptgetname(init) ) ;
  
  /* Build total transformation */
  gptconcattransform(&MN,&init->e,&N) ;
  init->e = MN ;

  gptaddEBelement( init, linecurrent_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}


static int linecurrent_sim(gptpar *par,double t,struct linecurrent_info *info)
{
  double x2py2, zpa, zma, zpa2, zma2, rp, rm, xyfac ;

  x2py2=X*X+Y*Y ;
  zpa=Z+info->L ;
  zma=Z-info->L ;
  zpa2=zpa*zpa ;
  zma2=zma*zma ;
  rp=sqrt(x2py2+zpa2) ;
  rm=sqrt(x2py2+zma2) ;
  if(rp==0 || rm==0) return( 1 ) ;
  xyfac=1.0/(rp*(zpa+rp)) - 1.0/(rm*(zma+rm)) ;

  BX = info->mfac*xyfac*Y ;
  BY =-info->mfac*xyfac*X ;
  BZ = 0.0 ;

  return( 1 ) ;
}

static int buildxyfromz(gpttransform *t)
{
  double x[3], y[3], z[3], temp[3] ;
  double xlen, ylen, zlen ;
  double zmin ;
  int i, imin ;

  /* Normalize z */
  for( i=0 ; i<3 ; i++ ) z[i] = t->m[i][2] ;
  zlen = gptVECLEN( z ) ;
  if( zlen==0 ) return( 1 ) ;
  for( i=0 ; i<3 ; i++ ) z[i] /= zlen ;

  /* Create temp vector */
  zmin = DBL_MAX ;
  imin = 0 ;
  for( i=0 ; i<3 ; i++ ) temp[i] = 0 ;
  for( i=0 ; i<3 ; i++ ) if( fabs(z[i])<zmin ) { zmin=fabs(z[i]) ; imin=i ; } ;
  temp[imin]=1 ;

  /* Calculate x vector */
  x[0] = z[1]*temp[2] - z[2]*temp[1] ;
  x[1] = z[2]*temp[0] - z[0]*temp[2] ;
  x[2] = z[0]*temp[1] - z[1]*temp[0] ;
  xlen = gptVECLEN( x ) ;
  for( i=0 ; i<3 ; i++ ) x[i] /= xlen ;

  /* Calculate y vector */
  y[0] = z[1]*x[2] - z[2]*x[1] ;
  y[1] = z[2]*x[0] - z[0]*x[2] ;
  y[2] = z[0]*x[1] - z[1]*x[0] ;
  ylen = gptVECLEN( y ) ;
  for( i=0 ; i<3 ; i++ ) y[i] /= ylen ;

  /* Build trasformation matrix */
  for( i=0 ; i<3 ; i++ ) t->m[i][0] = x[i] ;
  for( i=0 ; i<3 ; i++ ) t->m[i][1] = y[i] ;
  for( i=0 ; i<3 ; i++ ) t->m[i][2] = z[i] ;

  return( 0 ) ;
}
