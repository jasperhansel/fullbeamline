/* stdx.c - calculate std(x) */

#include "gdfa.h"

int stdx_func( double *result )
{
  /* Declarations */
  int num, tmpnum ;
  double *x, *nmacro ;

  /* Get selected arrays from GDFA kernel */
  if( gdfmgetarr( "x", &x, &num ) || num<1  ||
      gdfmgetarr( "nmacro", &nmacro, &tmpnum ) || tmpnum!=num ) return(1) ;

  /* Store result and return without error code */
  *result = gdfastd(nmacro,x,num) ;
  return( 0 ) ;
}
