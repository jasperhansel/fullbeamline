/* stdy.c - calculate std(y) */

#include "gdfa.h"

int stdy_func( double *result )
{
  /* Declarations */
  int num, tmpnum ;
  double *y, *nmacro ;

  /* Get selected arrays from GDFA kernel */
  if( gdfmgetarr( "y", &y, &num ) || num<1  ||
      gdfmgetarr( "nmacro", &nmacro, &tmpnum ) || tmpnum!=num ) return(1) ;

  /* Store result and return without error code */
  *result = gdfastd(nmacro,y,num) ;
  return( 0 ) ;
}
