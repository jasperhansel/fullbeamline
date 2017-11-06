/* TErectcavity.c - Rectangular TEmnp mode resonant cavity */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

struct tecav_info
{
  double a, b, d ;
  double ax, by, dz ;
  double fac, facoverkc2 ;
  double phi, w ;
} ;

static int tecav_sim(gptpar *par, double t, struct tecav_info *info) ;


void TErectcavity_init(gptinit *init)
{
  struct tecav_info *info ;
  int m, n, p ;
  double ax, by, dz ;
  double kc2 ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=8 )
    gpterror( "Syntax: %s(ECS,a,b,d,m,n,p,const,phi)\n", gptgetname(init) ) ;

  info = (struct tecav_info *)gptmalloc( sizeof(struct tecav_info) ) ;

  info->a = gptgetargdouble(init,1)/2 ;
  info->b = gptgetargdouble(init,2)/2 ;
  info->d = gptgetargdouble(init,3)/2 ;
  m = gptgetargint(init,4) ;
  n = gptgetargint(init,5) ;
  p = gptgetargint(init,6) ;
  info->fac = gptgetargdouble(init,7) ;
  info->phi = gptgetargdouble(init,8) ;

  if( info->a<=0 || info->b<=0 || info->d<=0 )
    gpterror( "Cavity dimensions must be larger than zero\n" ) ;

  if( m<0 || n<0 || p<0 )
    gpterror( "Eigenmode numbers must be positive or zero\n" ) ;

  if( p==0 )
    gpterror( "p eigenmode number may not be zero\n" ) ;

  if( m==0 && n==0 )
    gpterror( "m and n eigenmode numbers cannot both be zero\n" ) ;

  if( info->fac==0 )
    gptwarning( "Const is equal to zero, no power in cavity\n" ) ;

  info->ax  = ax = m*gpt_pi/(2*info->a) ;
  info->by  = by = n*gpt_pi/(2*info->b) ;
  info->dz  = dz = p*gpt_pi/(2*info->d) ;
  info->w   = sqrt(ax*ax+by*by+dz*dz)*gpt_c ;
  kc2 = ax*ax+by*by ;
  info->facoverkc2 = info->fac/kc2 ;  

  gptaddEBelement( init, tecav_sim, gptfree, GPTELEM_LOCAL, info ) ;
}


static int tecav_sim(gptpar *par, double t, struct tecav_info *info)
{
  double x, y, z ;  
  double ax, by, dz, w, facoverkc2 ;
  double sinax, cosax, sinby, cosby, sindz, cosdz, sinwt, coswt ;

  if( fabs(Z) > info->d ) return( 0 ) ;
  if( fabs(X) > info->a || fabs(Y) > info->b ) return( 0 ) ;

  x = X + info->a ;
  y = Y + info->b ;
  z = Z + info->d ;

  ax = info->ax ;
  by = info->by ;
  dz = info->dz ;
  w  = info->w  ;
  facoverkc2 = info->facoverkc2 ;

  sinax = sin(ax*x) ;
  cosax = cos(ax*x) ;
  sinby = sin(by*y) ;
  cosby = cos(by*y) ;
  sindz = sin(dz*z) ;
  cosdz = cos(dz*z) ;
  sinwt = sin(w*t+info->phi) ;
  coswt = cos(w*t+info->phi) ;

  BX = -facoverkc2*dz*ax*sinax*cosby*cosdz*sinwt ;
  BY = -facoverkc2*dz*by*cosax*sinby*cosdz*sinwt ;
  BZ =  info->fac*cosax*cosby*sindz*sinwt ;
  EX =  w*facoverkc2*by*cosax*sinby*sindz*coswt ;
  EY = -w*facoverkc2*ax*sinax*cosby*sindz*coswt ;

  return( 1 ) ;
}
