/* addzdiv.c - Add z-divergence (in eV/m) */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

void addzdiv_init(gptinit *init)
{
  gptparset *set ;
  gptinitpar *par ;
  char *name ;
  int i, len ;
  double div, zc, m0, Gold, Gnew, GBz2 ;

  if( gptgetargnum(init)!=3 )
    gpterror( "Syntax: %s(set,zc,div)\n", gptgetname(init) ) ;

  name = gptgetargstring(init,1) ;
  zc   = gptgetargdouble(init,2) ;
  div  = gptgetargdouble(init,3) ;

  /* Get particle set */
  if( gpttestparset( name )==NULL )
    gptwarning( "The particle set \"%s\" does not exist\n", name ) ;
  set = gptgetparset( name ) ;
  par = gptgetparsetpars( set,&len ) ;
  if( len==0 ) return ;

  /* Check particles */
  m0 = par[0].m ;
  for( i=0 ; i<len ; i++ ) if( par[i].m!=m0 )
   gpterror( "All particles must have identical mass\n" ) ;

  /* Add z-divergence */
  for( i=0 ; i<len ; i++ )
  {
    Gold = sqrt(1+gptVECSQR(par[i].GBr)) ;
    Gnew = Gold + div * (par[i].Wr[2]-zc) / (m0*gpt_c*gpt_c / -gpt_qe) ;
    if( Gnew < 1 )
      gpterror( "Adding divergence results in particle gamma<1\n" ) ;

    GBz2 = Gnew*Gnew - par[i].GBr[0]*par[i].GBr[0] - par[i].GBr[1]*par[i].GBr[1] - 1 ;
    if( GBz2 < 0 )
      gpterror( "Not able to maintain x- and y-momenta\n" ) ;

    par[i].GBr[2] = sqrt( GBz2 ) ;
  }
}
