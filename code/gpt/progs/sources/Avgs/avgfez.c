/* avgfEz.c - calculate average of fEz */

#include "gdfa.h"

int avgfEz_func( double *result )
{
  /* Declarations */
  int num, tmpnum ;
  double *fEz, *nmacro ;

  /* Get selected arrays from GDFA kernel */
  if( gdfmgetarr( "fEz", &fEz, &num ) || num<1  ||
      gdfmgetarr( "nmacro", &nmacro, &tmpnum ) || tmpnum!=num ) return(1) ;

  /* Store result and return without error code */
  *result = gdfamean(nmacro,fEz,num) ;
  return( 0 ) ;
}
