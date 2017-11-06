/* setGdist.c - Set particle energy distribution */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "elem.h"

void setGdist_init(gptinit *init)
{
  gptparset *set ;
  double *dist ;
  gptinitpar *par ;
  char *name ;
  int i, len ;
  double GBz2 ;

  if( gptgetargnum(init)<2 )
    gpterror( "Syntax: %s(set,DIST)\n", gptgetname(init) ) ;

  name = gptgetargstring(init,1) ;

  /* Get particle set */
  if( gpttestparset( name )==NULL )
    gptwarning( "The particle set \"%s\" does not exist\n", name ) ;
  set = gptgetparset( name ) ;
  par = gptgetparsetpars( set,&len ) ;

  /* Set distribution */
  dist = gptgetdistribution(init,len,2,gptdist_1D,set) ;
  for( i=0 ; i<len ; i++ )
  {
    GBz2 = dist[i]*dist[i] - par[i].GBr[0]*par[i].GBr[0] - par[i].GBr[1]*par[i].GBr[1] - 1 ;
    if( GBz2 < 0 )
      gpterror( "Not able to set gamma: Gamma^2 Betaz^2 becomes negative.\n" ) ;

    par[i].GBr[2] = sqrt( GBz2 ) ;
  }
  gptfree( dist ) ;
}
