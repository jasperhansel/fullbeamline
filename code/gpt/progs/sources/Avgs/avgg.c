/* avgG.c - calculate average G */

#include "gdfa.h"

int avgG_func( double *result )
{
  /* Declarations */
  int num, tmpnum ;
  double *G, *nmacro ;

  /* Get selected arrays from GDFA kernel */
  if( gdfmgetarr( "G", &G, &num ) || num<1  ||
      gdfmgetarr( "nmacro", &nmacro, &tmpnum ) || tmpnum!=num ) return(1) ;

  /* Store result and return without error code */
  *result = gdfamean(nmacro,G,num) ;
  return( 0 ) ;
}
