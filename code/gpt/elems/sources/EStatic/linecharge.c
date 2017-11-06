/* linecharge.c - Line with homogeneous charge density */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

#define ffac (1/(4*gpt_pi*gpt_eps0))

struct linecharge_info
{
  double mfac ;
  double L ;
} ;

static int linecharge_sim(gptpar *par,double t,struct linecharge_info *info) ;


void linecharge_init(gptinit *init)
{
  struct linecharge_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=2 )
    gpterror( "Syntax: %s(ECS,L,lambda)\n", gptgetname(init) ) ;

  info = (struct linecharge_info *)gptmalloc( sizeof(struct linecharge_info) ) ;

  info->L    = gptgetargdouble(init,1)/2 ;
  info->mfac = ffac*gptgetargdouble(init,2) ;

  gptaddEBelement( init, linecharge_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}


static int linecharge_sim(gptpar *par,double t,struct linecharge_info *info)
{
  double x2py2=X*X+Y*Y ;
  double zpa=Z+info->L ;
  double zma=Z-info->L ;
  double zpa2=zpa*zpa ;
  double zma2=zma*zma ;
  double rp=sqrt(x2py2+zpa2) ;
  double rm=sqrt(x2py2+zma2) ;
  double xyfac=1.0/(rm*(rm+zma))-1.0/(rp*(rp+zpa)) ;

  EX=info->mfac*xyfac*X ;
  EY=info->mfac*xyfac*Y ;
  EZ=info->mfac*(1.0/rm-1.0/rp) ;

  return( 1 ) ;
}

