/* avgt.c - calculate average time */

#include "gdfa.h"

int avgt_func( double *result )
{
  /* Declarations */
  int num, tmpnum ;
  double *t, *nmacro ;

  /* Get selected arrays from GDFA kernel */
  if( gdfmgetarr( "t", &t, &num ) || num<1  ||
      gdfmgetarr( "nmacro", &nmacro, &tmpnum ) || tmpnum!=num ) return(1) ;

  /* Store result and return without error code */
  *result = gdfamean(nmacro,t,num) ;
  return( 0 ) ;
}
