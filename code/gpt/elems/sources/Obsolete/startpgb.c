/* startpgb.c - Start a single particle using the "impuls" */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"


void startpgb_init(gptinit *init)
{
  double r[3], GBr[3] ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=6 )
    gpterror( "Syntax: %s(ECS,x,y,z,GBx,GBy,GBz)\n", gptgetname(init) ) ;

    r[0] = gptgetargdouble(init,1) ;
    r[1] = gptgetargdouble(init,2) ;
    r[2] = gptgetargdouble(init,3) ;
  GBr[0] = gptgetargdouble(init,4) ;
  GBr[1] = gptgetargdouble(init,5) ;
  GBr[2] = gptgetargdouble(init,6) ;

  gptaddparticleGB(init,r,GBr) ;
}
