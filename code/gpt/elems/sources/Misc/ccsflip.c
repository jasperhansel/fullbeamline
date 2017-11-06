/* ccsflip.c - Transfer a particle to an other axis */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include "elem.h"

struct ccsflip_info
{
  struct axis *toaxis ;
} ;

static int ccsflip_sim(gptpar *par,double t,struct ccsflip_info *info) ;


void ccsflip_init(gptinit *init)
{
  struct ccsflip_info *info ;
  char *toaxisname ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=1 )
    gpterror( "Syntax: %s(ECS,toaxis)\n", gptgetname(init) ) ;

  info = (struct ccsflip_info *)gptmalloc( sizeof(struct ccsflip_info) ) ;

  toaxisname = gptgetargstring(init,1) ;

  if( (info->toaxis=getaxis(toaxisname))==NULL )
    gpterror( "Specified axis \"%s\" not found\n", toaxisname ) ;

  gptaddEBelement( init, ccsflip_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}


static int ccsflip_sim(gptpar *par,double t,struct ccsflip_info *info)
{
  if( Z>0 ) par->newaxis = info->toaxis ;

  return( 1 ) ;
}
