/* avgfEy.c - calculate average of fEy */

#include "gdfa.h"

int avgfEy_func( double *result )
{
  /* Declarations */
  int num, tmpnum ;
  double *fEy, *nmacro ;

  /* Get selected arrays from GDFA kernel */
  if( gdfmgetarr( "fEy", &fEy, &num ) || num<1  ||
      gdfmgetarr( "nmacro", &nmacro, &tmpnum ) || tmpnum!=num ) return(1) ;

  /* Store result and return without error code */
  *result = gdfamean(nmacro,fEy,num) ;
  return( 0 ) ;
}
