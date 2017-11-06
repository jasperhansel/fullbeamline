/* dtmaxt.c - Enfore a dtmax within a specified time interval */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "elem.h"

struct dtmaxt_info
{
  double tstart ;
  double tend ;
  double dtmax ;
} ;

static int dtmaxt_out(double t, double *dt, double *x, void *vinfo) ;

void dtmaxt_init(gptinit *init)
{
  struct dtmaxt_info *info ;

  if( gptgetargnum(init)<3 )
    gpterror( "Syntax: %s(tstart,tend,dtmax)\n", gptgetname(init) ) ;

  info = (struct dtmaxt_info *)gptmalloc( sizeof(struct dtmaxt_info) ) ;

  info->tstart = gptgetargdouble(init,1) ;
  info->tend   = gptgetargdouble(init,2) ;
  info->dtmax  = gptgetargdouble(init,3) ;

  odeaddoutfunction( ODEFNC_USR, dtmaxt_out, info ) ;
  gptaddmainfunction( GPTMAINFNC_TER, gptfree, info ) ;
}


static int dtmaxt_out(double t, double *dt, double *x, void *vinfo)
{
  struct dtmaxt_info *info = (struct dtmaxt_info *)vinfo ;

  if( t>=info->tend ) return( 0 ) ;

  if( t>=info->tstart  ) {
    if( *dt>info->dtmax ) *dt=info->dtmax ;
  } else {
    if( t+*dt>info->tstart+info->dtmax ) *dt=info->tstart+info->dtmax - t ;
  }

  return( 0 ) ;
}
