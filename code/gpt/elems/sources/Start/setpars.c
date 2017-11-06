/* setpars.c - Set number, mass and charge of particles in a new set */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"


void setparticles_init(gptinit *init)
{
  gptparset *set ;
  char *name ;
  double m, q, Qtot ;
  int i, N ;

  double Wr[3]   = {0.0, 0.0, 0.0} ;
  double WGBr[3] = {0.0, 0.0, 0.0} ;

  if( gptgetargnum(init)!=5 )
    gpterror( "Syntax: %s(set,N,m,q,Qtot)\n", gptgetname(init) ) ;

  name = gptgetargstring(init,1) ;
  N    = gptgetargint   (init,2) ;
  m    = gptgetargdouble(init,3) ;
  q    = gptgetargdouble(init,4) ;
  Qtot = gptgetargdouble(init,5) ;

  if( N<0 )
    gpterror( "Invalid number of particles\n" ) ;
  if( N==0 )
    gptwarning( "Empty particle set\n" ) ;
  if( gpttestparset( name )!=NULL )
    gpterror( "The particle set %s already exists\n", name ) ;
  if( q*Qtot<0.0 )
    gpterror( "q and Qtot have opposite sign\n" ) ;

  set = gptgetparset( name ) ;

  for(i=0 ; i<N ; i++) gptaddparmqn( set, Wr, WGBr, m, q, Qtot/q/N ) ;
}
