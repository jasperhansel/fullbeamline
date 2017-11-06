/* setstrli.c - Start a number of particles on straight line */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

void setstartline_init(gptinit *init)
{
  int i,j, numarg, N ;
  double r[3], r1[3], r2[3], GBr[3] ;
  double m,q,n ;
  gptparset *set ;

  numarg = gptgetargnum(init) ;
  if( numarg!=8 && numarg!=11 && numarg!=14)
    gpterror( "Syntax: %s(set,N,x1,y1,z1,x2,y2,z2,[GBx,GBy,GBz,[m,q,n]])\n", gptgetname(init) ) ;

     set = gptgetparset( gptgetargstring(init,1) ) ;
       N = gptgetargint(init,2) ;
   r1[0] = gptgetargdouble(init,3) ;
   r1[1] = gptgetargdouble(init,4) ;
   r1[2] = gptgetargdouble(init,5) ;
   r2[0] = gptgetargdouble(init,6) ;
   r2[1] = gptgetargdouble(init,7) ;
   r2[2] = gptgetargdouble(init,8) ;

  if( N<2 )
    gpterror( "%s: N must be at least 2\n", gptgetname(init) ) ;

  if( numarg> 8 )
  {
    GBr[0] = gptgetargdouble(init, 9) ;
    GBr[1] = gptgetargdouble(init,10) ;
    GBr[2] = gptgetargdouble(init,11) ;
  } else
    GBr[0] = GBr[1] = GBr[2] = 0.0 ;

  if( numarg<=11 )
  {
    for(i=0 ; i<N ; i++) 
    {
      for(j=0 ; j<3 ; j++) r[j] = r1[j] + (double)i*(r2[j]-r1[j])/(N-1) ;
      gptaddpar(set,r,GBr) ;
    }
  } else
  {
    m = gptgetargdouble(init,12) ;
    q = gptgetargdouble(init,13) ;
    n = gptgetargdouble(init,14) ;

    for(i=0 ; i<N ; i++) 
    {
      for(j=0 ; j<3 ; j++) r[j] = r1[j] + (double)i*(r2[j]-r1[j])/(N-1) ;
      gptaddparmqn(set,r,GBr,m,q,n) ;
    }
  }
}
