/* avgz.c - calculate average z */

#include "gdfa.h"

int avgz_func( double *result )
{
  /* Declarations */
  int num, tmpnum ;
  double *z, *nmacro ;

  /* Get selected arrays from GDFA kernel */
  if( gdfmgetarr( "z", &z, &num ) || num<1  ||
      gdfmgetarr( "nmacro", &nmacro, &tmpnum ) || tmpnum!=num ) return(1) ;

  /* Store result and return without error code */
  *result = gdfamean(nmacro,z,num) ;
  return( 0 ) ;
}
