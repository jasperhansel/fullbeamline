/* avgfBy.c - calculate average of fBy */

#include "gdfa.h"

int avgfBy_func( double *result )
{
  /* Declarations */
  int num, tmpnum ;
  double *fBy, *nmacro ;

  /* Get selected arrays from GDFA kernel */
  if( gdfmgetarr( "fBy", &fBy, &num ) || num<1  ||
      gdfmgetarr( "nmacro", &nmacro, &tmpnum ) || tmpnum!=num ) return(1) ;

  /* Store result and return without error code */
  *result = gdfamean(nmacro,fBy,num) ;
  return( 0 ) ;
}
