/* avgfBz.c - calculate average of fBz */

#include "gdfa.h"

int avgfBz_func( double *result )
{
  /* Declarations */
  int num, tmpnum ;
  double *fBz, *nmacro ;

  /* Get selected arrays from GDFA kernel */
  if( gdfmgetarr( "fBz", &fBz, &num ) || num<1  ||
      gdfmgetarr( "nmacro", &nmacro, &tmpnum ) || tmpnum!=num ) return(1) ;

  /* Store result and return without error code */
  *result = gdfamean(nmacro,fBz,num) ;
  return( 0 ) ;
}
