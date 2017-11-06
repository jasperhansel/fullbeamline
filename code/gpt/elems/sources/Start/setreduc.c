/* setreduc.c - Reduce the number of particles in a set */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"


void setreduce_init(gptinit *init)
{
  gptparset *set ;
  gptinitpar *par, tmppar ;
  char *name ;
  int i,j, len, N ;

  if( gptgetargnum(init)!=2 )
    gpterror( "Syntax: %s(set,N)\n", gptgetname(init) ) ;

  name = gptgetargstring(init,1) ;
  N    = gptgetargint   (init,2) ;

  /* Get particle set */
  if( gpttestparset( name )==NULL )
    gptwarning( "The particle set \"%s\" does not exist\n", name ) ;
  set = gptgetparset( name ) ;
  par = gptgetparsetpars( set,&len ) ;

  /* Check number of particles */
  if( N>len )
    gpterror( "Reduction from %d to %d particles is not possible.\n", len, N ) ;

  /* Select first N particles randomly */
  for( i=0 ; i<N ; i++ )
  {
    j = (rand()%(len-i)) ;
    tmppar = par[len-j-1] ;
    par[len-j-1] = par[i] ;
    par[i] = tmppar ;
  }

  /* Adapt total charge and remove remaining particles*/
  for( i=0 ; i<N ; i++ ) par[i].n *= (double)len/(double)N ;
  while( len!=N ) gptremovepar(set,len-1,&len ) ;
}
