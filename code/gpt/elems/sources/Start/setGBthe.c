/* setGBthe.c - Set angle distribution in GBz-GBxy plane */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "elem.h"

void setGBthetadist_init(gptinit *init)
{
  gptparset *set ;
  double *dist ;
  gptinitpar *par ;
  char *name ;
  int i, len ;
  double GBrxyz ;

  if( gptgetargnum(init)<2 )
    gpterror( "Syntax: %s(set,DIST)\n", gptgetname(init) ) ;

  name = gptgetargstring(init,1) ;

  /* Get particle set */
  if( gpttestparset( name )==NULL )
    gptwarning( "The particle set \"%s\" does not exist\n", name ) ;
  set = gptgetparset( name ) ;
  par = gptgetparsetpars( set,&len ) ;

  /* Set distribution */
  dist = gptgetdistribution(init,len,2,gptdist_2DSin,set) ;
  for( i=0 ; i<len ; i++ )
  {
    if( par[i].GBr[0]==0.0 && par[i].GBr[1]==0.0 )
    {
      /* Make phi angle 0 */
      par[i].GBr[0] = par[i].GBr[2]*sin(dist[i]) ;
      par[i].GBr[2] *= cos(dist[i]) ;
    } else
    {
      /* Keep existing phi angle */
      GBrxyz = sqrt( par[i].GBr[0]*par[i].GBr[0] + par[i].GBr[1]*par[i].GBr[1] + par[i].GBr[2]*par[i].GBr[2] ) ;
      par[i].GBr[2]  = GBrxyz*cos(dist[i]) ;
      gptr2carth(GBrxyz*sin(dist[i]),par[i].GBr[0],par[i].GBr[1],&par[i].GBr[0],&par[i].GBr[1]) ;
    }
  }
  gptfree( dist ) ;
}
