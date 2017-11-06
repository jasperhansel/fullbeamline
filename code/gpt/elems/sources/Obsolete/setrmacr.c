/* setrmacr.c - Set the radius of all particles in a set */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"


void setrmacro_init(gptinit *init)
{
  gptparset *set ;
  gptinitpar *par ;
  char *name ;
  int i, len ;
  double R ;

  if( gptgetargnum(init)!=2 )
    gpterror( "Syntax: %s(set,R)\n", gptgetname(init) ) ;

  name = gptgetargstring(init,1) ;
  R    = gptgetargdouble(init,2) ;

  if( gpttestparset( name )==NULL )
    gptwarning( "The particle set %s does not exist\n", name ) ;

  if( R<=0.0 )
    gpterror( "The specified radius of the particles (%g) must be positive\n", R ) ;

  set = gptgetparset( name ) ;
  par = gptgetparsetpars( set,&len ) ;

  for( i=0 ; i<len ; i++ ) par[i].r = R ;
}
