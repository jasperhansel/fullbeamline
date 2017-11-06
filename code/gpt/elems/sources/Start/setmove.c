/* setmove.c - Move all particles within a set */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"


void setmove_init(gptinit *init)
{
  gptparset *setfrom, *setto ;
  gptinitpar *pars ;
  char *namefrom, *nameto ;
  int i, lenfrom ;

  if( gptgetargnum(init)!=2 )
    gpterror( "Syntax: %s(setfrom,setto)\n", gptgetname(init) ) ;

  namefrom = gptgetargstring(init,1) ;
  nameto   = gptgetargstring(init,2) ;

  /* Get particle set from */
  if( gpttestparset( namefrom )==NULL )
    gptwarning( "The particle set \"%s\" does not exist\n", namefrom ) ;
  setfrom = gptgetparset( namefrom ) ;
  pars = gptgetparsetpars( setfrom, &lenfrom ) ;

  /* Get particle set to */
  setto = gptgetparset( nameto ) ;

  /* Copy from to to */
  for( i=0 ; i<lenfrom ; i++ )
	  gptaddparmqnartid(setto,pars[i].Wr,pars[i].GBr,pars[i].m,pars[i].q,pars[i].n,pars[i].paxis,pars[i].r,pars[i].tstart,pars[i].ID) ;

  /* Set Hammersley dimension number to largest of from and to */
  if( setfrom->ndisthammersley > setto->ndisthammersley )
    setto->ndisthammersley = setfrom->ndisthammersley ;

  /* Remove original set */
  gptremoveparset(setfrom) ;
}
