/* CSalphax.c - Calculate Courant Snyder alpha for x x' subspace */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "gdfa.h"

int CSalphax_func( double *result )
{
  int num, tmpnum ;
  double *nmacro, *x,  *Bx, *Bz, *G  ;
  double *xc, *xpc ;
  double xx, xpxp, xxp, emix ;

  if( gdfmgetarr( "x" , &x,  &num    ) ||  num <=2     ||
      gdfmgetarr( "nmacro", &nmacro, &tmpnum ) || tmpnum!=num ||
      gdfmgetarr( "Bx", &Bx, &tmpnum ) ||  tmpnum!=num ||
      gdfmgetarr( "Bz", &Bz, &tmpnum ) ||  tmpnum!=num ||
      gdfmgetarr( "G" , &G,  &tmpnum ) ||  tmpnum!=num ) return(1) ;

  /* Allocate memory for new items xc, xpc */
  xc = (double *)malloc( num*sizeof(double)) ;
  xpc= (double *)malloc( num*sizeof(double)) ;
  if (xc==NULL || xpc==NULL ) { fprintf( stderr, "Not enough memory\n") ; return (1) ; }

  gdfasubavg(nmacro,xc, x, num) ;
  gdfasubavg(nmacro,xpc,Bx,num) ;

  xx   = gdfamean2(nmacro,xc ,xc ,num) ;
  xpxp = gdfamean2(nmacro,xpc,xpc,num) ;
  xxp  = gdfamean2(nmacro,xc ,xpc,num) ;

  emix = stdsqrt( xx*xpxp-xxp*xxp) ;
  if( emix==0.0 ) return( 1 ) ;

  *result = -xxp/emix ;

  free( xc ) ;
  free( xpc ) ;

  return(0) ;
}
