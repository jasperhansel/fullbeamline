/* viscosity.c:  */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

struct viscosity_info
{
  double eta ;
} ;

static void viscosity_sim( double t, double *x, double *p, void *vinfo ) ;

void viscosity_init(gptinit *init)
{
  struct viscosity_info *info ;

  if( gptgetargnum(init)!=1 )
    gpterror( "Syntax: %s(eta)\n", gptgetname(init) ) ;

  info = (struct viscosity_info *)gptmalloc( sizeof(struct viscosity_info) ) ;

  info->eta = gptgetargdouble(init,1) ;

  /* Register viscosity function */
  odeaddfprfunction( ODEFNC_INT, viscosity_sim,info ) ;
}


static void viscosity_sim( double t, double *x, double *p, void *vinfo )
{
  struct par *par ;
  double Br[3], fac ;
  struct viscosity_info *info = (struct viscosity_info *)vinfo ;

  for(int i=0 ; i<numpar ; i++)
  {
    par = &pars[i] ;

    Br[0] = par->GBr[0]/par->G ;
    Br[1] = par->GBr[1]/par->G ;
    Br[2] = par->GBr[2]/par->G ;

    fac = -6*gpt_pi*info->eta*sqrt(par->r2)*gpt_c ;

    par->WF[0] += fac*Br[0] ;
    par->WF[1] += fac*Br[1] ;
    par->WF[2] += fac*Br[2] ;
  }
}
