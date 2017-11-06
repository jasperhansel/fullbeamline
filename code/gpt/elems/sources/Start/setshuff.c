/* setshuffle.c - Shuffle all particles within a set */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"


void setshuffle_init(gptinit *init)
{
  gptparset *set ;
  gptinitpar *par, tmppar ;
  char *name ;
  int i,j, len ;

  if( gptgetargnum(init)!=1 )
    gpterror( "Syntax: %s(set)\n", gptgetname(init) ) ;

  name = gptgetargstring(init,1) ;

  /* Get particle set */
  if( gpttestparset( name )==NULL )
    gptwarning( "The particle set \"%s\" does not exist\n", name ) ;
  set = gptgetparset( name ) ;
  par = gptgetparsetpars( set,&len ) ;

  /* Shuffle */
  for( i=0 ; i<len ; i++ )
  {
    j = (intpprand()%(len-i)) ;
    tmppar = par[len-j-1] ;
    par[len-j-1] = par[i] ;
    par[i] = tmppar ;
  }
}
