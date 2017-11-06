/* addydiv.c - Add y-divergence (in rad/meter) */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"


void addydiv_init(gptinit *init)
{
  gptparset *set ;
  gptinitpar *par ;
  char *name ;
  int i, len ;
  double yc, div, Gold, GBz2 ;

  if( gptgetargnum(init)!=3 )
    gpterror( "Syntax: %s(set,yc,div)\n", gptgetname(init) ) ;

  name = gptgetargstring(init,1) ;
  yc   = gptgetargdouble(init,2) ;
  div  = gptgetargdouble(init,3) ;

  /* Get particle set */
  if( gpttestparset( name )==NULL )
    gptwarning( "The particle set \"%s\" does not exist\n", name ) ;
  set = gptgetparset( name ) ;
  par = gptgetparsetpars( set,&len ) ;

  /* Add y-divergence */
  for( i=0 ; i<len ; i++ )
  {
    Gold = sqrt(1+gptVECSQR(par[i].GBr)) ;

    par[i].GBr[1] += div*(par[i].Wr[1]-yc) ;

    GBz2 = Gold*Gold - par[i].GBr[0]*par[i].GBr[0] - par[i].GBr[1]*par[i].GBr[1] - 1 ;
    if( GBz2 < 0 )
      gpterror( "Not able to maintain gamma\n" ) ;

    par[i].GBr[2] = sqrt( GBz2 ) ;
  }
}
