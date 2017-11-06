/* stdBz.c - calculate std(Bz) */

#include "gdfa.h"

int stdBz_func( double *result )
{
  /* Declarations */
  int num, tmpnum ;
  double *Bz, *nmacro ;

  /* Get selected arrays from GDFA kernel */
  if( gdfmgetarr( "Bz", &Bz, &num ) || num<1  ||
      gdfmgetarr( "nmacro", &nmacro, &tmpnum ) || tmpnum!=num ) return(1) ;

  /* Store result and return without error code */
  *result = gdfastd(nmacro,Bz,num) ;
  return( 0 ) ;
}
