/* stdBy.c - calculate std(By) */

#include "gdfa.h"

int stdBy_func( double *result )
{
  /* Declarations */
  int num, tmpnum ;
  double *By, *nmacro ;

  /* Get selected arrays from GDFA kernel */
  if( gdfmgetarr( "By", &By, &num ) || num<1  ||
      gdfmgetarr( "nmacro", &nmacro, &tmpnum ) || tmpnum!=num ) return(1) ;

  /* Store result and return without error code */
  *result = gdfastd(nmacro,By,num) ;
  return( 0 ) ;
}
