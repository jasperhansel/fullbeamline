/* cputime.c - get cputime */

#include "gdfa.h"

int cputime_func( double *result ) 
{
  if( gdfmgetval( "cputime", result ) ) return( 1 ) ;
  return( 0 ) ;
}
