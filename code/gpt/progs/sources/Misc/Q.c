/* Q.c: Total charge */

#include "gdfa.h"

int Q_func( double *result )
{
  /* Declarations */
  int i, num, tmpnum ;
  double *nmacro, *q ;
  double sumnq ;

  /* Get selected arrays from GDFA kernel */
  if( gdfmgetarr( "nmacro", &nmacro, &num ) || num<1  ||
      gdfmgetarr( "q", &q, &tmpnum ) || tmpnum!=num ) return(1) ;

  /* Calculate total charge present in simulation */
  sumnq = 0.0 ;
  for(i=0 ; i<num ; i++) sumnq += nmacro[i]*q[i] ;

  /* Store result and return without error code */
  *result = sumnq ;
  return( 0 ) ;
}
