/* magplate.c - Rectangular plate with homogeneous magnetic charge density */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

#define ffac (gpt_mu0/(4*gpt_pi))

struct magplate_info
{
  double mfac ;
  double a ;
  double b ;
} ;

static int magplate_sim(gptpar *par,double t,struct magplate_info *info) ;


void magplate_init(gptinit *init)
{
  struct magplate_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=3 )
    gpterror( "Syntax: %s(ECS,a,b,sigma)\n", gptgetname(init) ) ;

  info = (struct magplate_info *)gptmalloc( sizeof(struct magplate_info) ) ;

  info->a    = gptgetargdouble(init,1)/2 ;
  info->b    = gptgetargdouble(init,2)/2 ;
  info->mfac = ffac*gptgetargdouble(init,3) ;

  gptaddEBelement( init, magplate_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}


static void tBplate(double x, double y, double z, int sign, double *B)
{
  double r=sqrt(x*x+y*y+z*z) ;
  B[0]-=sign*log(r+y) ;
  B[1]-=sign*log(r+x) ;
  B[2]+=sign*atan((x*y)/(r*z)) ;
}


static int magplate_sim(gptpar *par,double t,struct magplate_info *info)
{
  tBplate(X-info->a,Y-info->b,Z,+1,par->B) ;
  tBplate(X-info->a,Y+info->b,Z,-1,par->B) ;
  tBplate(X+info->a,Y-info->b,Z,-1,par->B) ;
  tBplate(X+info->a,Y+info->b,Z,+1,par->B) ;

  BX*=info->mfac ;
  BY*=info->mfac ;
  BZ*=info->mfac ;

  return( 1 ) ;
}
