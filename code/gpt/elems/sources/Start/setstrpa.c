/* setstrpa.c - Start a single particle in a set */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

void setstartpar_init(gptinit *init)
{
  int numarg ;
  double r[3], GBr[3] ;
  double m,q,n ;
  gptparset *set ;

  numarg = gptgetargnum(init) ;
  if( numarg!=7 && numarg!=10 )
    gpterror( "Syntax: %s(set,x,y,z,GBx,GBy,GBz,[m,q,n])\n", gptgetname(init) ) ;

     set = gptgetparset( gptgetargstring(init,1) ) ;
    r[0] = gptgetargdouble(init,2) ;
    r[1] = gptgetargdouble(init,3) ;
    r[2] = gptgetargdouble(init,4) ;
  GBr[0] = gptgetargdouble(init,5) ;
  GBr[1] = gptgetargdouble(init,6) ;
  GBr[2] = gptgetargdouble(init,7) ;

  if( numarg==7 )
  {
    gptaddpar(set,r,GBr) ;
  } else
  {
    m = gptgetargdouble(init,8) ;
    q = gptgetargdouble(init,9) ;
    n = gptgetargdouble(init,10) ;

	if( n<0 )
		gptwarning("Macro particle represents a negative number of elementary particles") ;

    gptaddparmqn(set,r,GBr,m,q,n) ;
  }
}
