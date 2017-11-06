/* setphidi.c - Set angle in xy-plane*/

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "elem.h"

void setphidist_init(gptinit *init)
{
  gptparset *set ;
  double *dist ;
  gptinitpar *par ;
  char *name ;
  int i, len ;
  double rxy ;

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
    rxy = sqrt( par[i].Wr[0]*par[i].Wr[0] + par[i].Wr[1]*par[i].Wr[1] ) ;
    par[i].Wr[0] = rxy*cos(dist[i]) ;
    par[i].Wr[1] = rxy*sin(dist[i]) ;
  }
  gptfree( dist ) ;
}
