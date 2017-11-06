/* Rmax.c - calculate max radius */

#include "gdfa.h"
#include <math.h>

int rmax_func( double *result ) 
{
  double *x, *y, r2;
  int numx, numy, i;

  *result = 0.0 ;
  if( gdfmgetarr( "x", &x, &numx ) || numx<2 ) return( 1 ) ;
      gdfmgetarr( "y", &y, &numy ) ;
  for(i=0 ; i<numx ; i++) 
    { r2 = x[i]*x[i] + y[i]*y[i] ;
      if (*result < r2 ) *result = r2 ;
    } ;
  *result = sqrt(*result) ;
  return( 0 ) ;
}

