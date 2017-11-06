/* barmagnet.c - Rectangular magnet with homogeneous magnetization */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

struct barmagnet_info
{
  double a,b,L ;
  double mfac ;
} ;

static int barmagnet_sim(gptpar *par,double t,struct barmagnet_info *info) ;


void barmagnet_init(gptinit *init)
{
  struct barmagnet_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=4 )
    gpterror( "Syntax: %s(ECS,a,b,L,B)\n", gptgetname(init) ) ;

  info = (struct barmagnet_info *)gptmalloc( sizeof(struct barmagnet_info) ) ;

  info->a    = gptgetargdouble(init,1)/2 ;
  info->b    = gptgetargdouble(init,2)/2 ;
  info->L    = gptgetargdouble(init,3)/2 ;
  info->mfac = gptgetargdouble(init,4)/(4*gpt_pi) ;

  gptaddEBelement( init, barmagnet_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}


static void tBplate(double x, double y, double z, int sign, double *B)
{
  double r=sqrt(x*x+y*y+z*z) ;
  B[0]-=sign*log(r+y) ;
  B[1]-=sign*log(r+x) ;
  B[2]+=sign*atan((x*y)/(r*z)) ;
}


static int barmagnet_sim(gptpar *par,double t,struct barmagnet_info *info)
{
  tBplate(X-info->a,Y-info->b,Z-info->L,+1,par->B) ;
  tBplate(X-info->a,Y+info->b,Z-info->L,-1,par->B) ;
  tBplate(X+info->a,Y-info->b,Z-info->L,-1,par->B) ;
  tBplate(X+info->a,Y+info->b,Z-info->L,+1,par->B) ;

  tBplate(X-info->a,Y-info->b,Z+info->L,-1,par->B) ;
  tBplate(X-info->a,Y+info->b,Z+info->L,+1,par->B) ;
  tBplate(X+info->a,Y-info->b,Z+info->L,+1,par->B) ;
  tBplate(X+info->a,Y+info->b,Z+info->L,-1,par->B) ;

  if( X<info->a && X>-info->a &&
      Y<info->b && Y>-info->b &&
      Z<info->L && Z>-info->L )
    BZ += 4*gpt_pi ;

  BX*=info->mfac ;
  BY*=info->mfac ;
  BZ*=info->mfac ;

  return( 1 ) ;
}
