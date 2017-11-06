/* magline.c - Line with homogeneous magnetic charge density */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

#define ffac (gpt_mu0/(4*gpt_pi))

struct magline_info
{
  double mfac ;
  double L ;
} ;

static int magline_sim(gptpar *par,double t,struct magline_info *info) ;


void magline_init(gptinit *init)
{
  struct magline_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=2 )
    gpterror( "Syntax: %s(ECS,L,lambda)\n", gptgetname(init) ) ;

  info = (struct magline_info *)gptmalloc( sizeof(struct magline_info) ) ;

  info->L    = gptgetargdouble(init,1)/2 ;
  info->mfac = ffac*gptgetargdouble(init,2) ;

  gptaddEBelement( init, magline_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}


static int magline_sim(gptpar *par,double t,struct magline_info *info)
{
  double x2py2=X*X+Y*Y ;
  double zpa=Z+info->L ;
  double zma=Z-info->L ;
  double zpa2=zpa*zpa ;
  double zma2=zma*zma ;
  double rp=sqrt(x2py2+zpa2) ;
  double rm=sqrt(x2py2+zma2) ;
  double xyfac=1.0/(rm*(rm+zma))-1.0/(rp*(rp+zpa)) ;

  BX=info->mfac*xyfac*X ;
  BY=info->mfac*xyfac*Y ;
  BZ=info->mfac*(1.0/rm-1.0/rp) ;

  return( 1 ) ;
}

