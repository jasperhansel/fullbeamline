/* setGBphi.c - Set angle distribution in GBx-GBy plane */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "elem.h"

void setGBphidist_init(gptinit *init)
{
  gptparset *set ;
  double *dist ;
  gptinitpar *par ;
  char *name ;
  int i, len ;
  double rGBxy ;

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
    rGBxy = sqrt( par[i].GBr[0]*par[i].GBr[0] + par[i].GBr[1]*par[i].GBr[1] ) ;
    par[i].GBr[0] = rGBxy*cos(dist[i]) ;
    par[i].GBr[1] = rGBxy*sin(dist[i]) ;
  }
  gptfree( dist ) ;
}

