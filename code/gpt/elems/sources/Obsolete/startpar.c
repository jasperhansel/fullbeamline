/* startpar.c - Start a single particle */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"


void startpar_init(gptinit *init)
{
  double r[3], Br[3] ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=6 )
    gpterror( "Syntax: %s(ECS,x,y,z,Bx,By,Bz)\n", gptgetname(init) ) ;

   r[0] = gptgetargdouble(init,1) ;
   r[1] = gptgetargdouble(init,2) ;
   r[2] = gptgetargdouble(init,3) ;
  Br[0] = gptgetargdouble(init,4) ;
  Br[1] = gptgetargdouble(init,5) ;
  Br[2] = gptgetargdouble(init,6) ;

  gptaddparticle(init,r,Br) ;
}
