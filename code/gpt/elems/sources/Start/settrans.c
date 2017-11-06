/* settrans.c: Transform a particle set */

#include <stdio.h>
#include <math.h>
#include "elem.h"

void settransform_init(gptinit *init)
{
  gptparset *set ;
  gptinitpar *par ;
  char *name ;
  int i, len ;
  double tmp[3] ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=1 )
    gpterror( "Syntax: %s(ecs,set)\n", gptgetname(init) ) ;

  name = gptgetargstring(init,1) ;

  /* Get particle set */
  if( gpttestparset( name )==NULL )
    gptwarning( "The particle set \"%s\" does not exist\n", name ) ;
  set = gptgetparset( name ) ;
  par = gptgetparsetpars( set,&len ) ;

  /* Transform every particle to ecs */
  for(i=0 ; i<len ; i++)
  {
    /* Transform to ecs */
    gpttoWCS( &init->e, par[i].Wr, tmp ) ;
    gpttoWCS( &init->paxis->a, tmp, par[i].Wr ) ;
    gptdirectiontoWCS( &init->e, par[i].GBr, tmp ) ;
    gptdirectiontoWCS( &init->paxis->a, tmp, par[i].GBr ) ;

    /* Store axis */
    par[i].paxis = init->paxis ;
  }
}
