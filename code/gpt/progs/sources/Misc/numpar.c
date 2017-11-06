/* numpar.c - calculate number of particles */

#include "gdfa.h"

int numpar_func( double *result ) 
{
  double *x ;
  int numx ;

  if( gdfmgetarr( "x", &x, &numx ) ) *result=0.0 ; else *result=numx ;
  return( 0 ) ;
}
