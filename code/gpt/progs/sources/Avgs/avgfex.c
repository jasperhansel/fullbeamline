/* avgfEx.c - calculate average of fEx */

#include "gdfa.h"

int avgfEx_func( double *result )
{
  /* Declarations */
  int num, tmpnum ;
  double *fEx, *nmacro ;

  /* Get selected arrays from GDFA kernel */
  if( gdfmgetarr( "fEx", &fEx, &num ) || num<1  ||
      gdfmgetarr( "nmacro", &nmacro, &tmpnum ) || tmpnum!=num ) return(1) ;

  /* Store result and return without error code */
  *result = gdfamean(nmacro,fEx,num) ;
  return( 0 ) ;
}
