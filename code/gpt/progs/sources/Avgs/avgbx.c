/* avgBx.c - calculate average Bx */

#include "gdfa.h"

int avgBx_func( double *result )
{
  /* Declarations */
  int num, tmpnum ;
  double *Bx, *nmacro ;

  /* Get selected arrays from GDFA kernel */
  if( gdfmgetarr( "Bx", &Bx, &num ) || num<1  ||
      gdfmgetarr( "nmacro", &nmacro, &tmpnum ) || tmpnum!=num ) return(1) ;

  /* Store result and return without error code */
  *result = gdfamean(nmacro,Bx,num) ;
  return( 0 ) ;
}
