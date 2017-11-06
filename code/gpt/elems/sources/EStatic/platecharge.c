/* platecharge.c - Rectangular plate with homogeneous charge density */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

#define ffac (1/(4*gpt_pi*gpt_eps0))

struct platecharge_info
{
  double mfac ;
  double a ;
  double b ;
} ;

static int platecharge_sim(gptpar *par,double t,struct platecharge_info *info) ;


void platecharge_init(gptinit *init)
{
  struct platecharge_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=3 )
    gpterror( "Syntax: %s(ECS,a,b,sigma)\n", gptgetname(init) ) ;

  info = (struct platecharge_info *)gptmalloc( sizeof(struct platecharge_info) ) ;

  info->a    = gptgetargdouble(init,1)/2 ;
  info->b    = gptgetargdouble(init,2)/2 ;
  info->mfac = ffac*gptgetargdouble(init,3) ;

  gptaddEBelement( init, platecharge_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}


static void tEplate(double x, double y, double z, int sign, double *E)
{
  double r=sqrt(x*x+y*y+z*z) ;
  E[0]-=sign*log(r+y) ;
  E[1]-=sign*log(r+x) ;
  E[2]+=sign*atan((x*y)/(r*z)) ;
}


static int platecharge_sim(gptpar *par,double t,struct platecharge_info *info)
{
  tEplate(X-info->a,Y-info->b,Z,+1,par->E) ;
  tEplate(X-info->a,Y+info->b,Z,-1,par->E) ;
  tEplate(X+info->a,Y-info->b,Z,-1,par->E) ;
  tEplate(X+info->a,Y+info->b,Z,+1,par->E) ;

  EX*=info->mfac ;
  EY*=info->mfac ;
  EZ*=info->mfac ;

  return( 1 ) ;
}


