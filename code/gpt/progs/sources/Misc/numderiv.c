/* numderivs.c - get number of calculated derivatives */

#include "gdfa.h"

int numderivs_func( double *result ) 
{
  if( gdfmgetval( "numderivs", result ) ) return( 1 ) ;
  return( 0 ) ;
}
