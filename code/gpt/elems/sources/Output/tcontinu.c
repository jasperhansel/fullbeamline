/* tcontinu.c: Continue the simulation till the specified time */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

struct tcontinue_info
{
  double t ;
} ;

static int tcontinue_out(double t, double *dt, double *x, void *vinfo) ;

void tcontinue_init(gptinit *init)
{
  struct tcontinue_info *info ;

  if( gptgetargnum(init)!=1 )
    gpterror( "Syntax: %s(t)\n", gptgetname(init) ) ;

  info = (struct tcontinue_info *)gptmalloc( sizeof(struct tcontinue_info) ) ;
  info->t        = gptgetargdouble(init,1) ;

  odeaddoutfunction( ODEFNC_USR, tcontinue_out, info ) ;
}

static int tcontinue_out(double t, double *dt, double *x, void *vinfo)
{
  struct tcontinue_info *info = (struct tcontinue_info *)vinfo ;

  if( t<info->t ) return(1) ;
  return(0) ;
}
