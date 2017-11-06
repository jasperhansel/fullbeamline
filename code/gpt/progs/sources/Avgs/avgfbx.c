/* avgfBx.c - calculate average of fBx */

#include "gdfa.h"

int avgfBx_func( double *result )
{
  /* Declarations */
  int num, tmpnum ;
  double *fBx, *nmacro ;

  /* Get selected arrays from GDFA kernel */
  if( gdfmgetarr( "fBx", &fBx, &num ) || num<1  ||
      gdfmgetarr( "nmacro", &nmacro, &tmpnum ) || tmpnum!=num ) return(1) ;

  /* Store result and return without error code */
  *result = gdfamean(nmacro,fBx,num) ;
  return( 0 ) ;
}
