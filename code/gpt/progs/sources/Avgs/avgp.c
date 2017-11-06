/* avgp.c - Calculate mean momentum of beam in [eV s / m] */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "gdfa.h"

int avgp_func( double *result )
{
  /* Declarations */
  int i, num, tmpnum ;
  double *GB, *G, *m, *nmacro ;
  double Beta0 ;

  /* Get selected arrays from GDFA kernel */
  if( gdfmgetarr( "G", &G, &num ) || num<1  ||
      gdfmgetarr( "m", &m, &tmpnum ) || tmpnum!=num ||
      gdfmgetarr( "nmacro", &nmacro, &tmpnum ) || tmpnum!=num ) return(1) ;

  /* Calculate Beta0=avg(GB)/avg(G) */
  GB = (double *)malloc( num*sizeof(double)) ;
  for(i=0 ; i<num ; i++) GB[i] = stdsqrt(G[i]*G[i]-1) ;
  Beta0 = gdfamean(nmacro,GB,num) / gdfamean(nmacro,G,num) ;
  free( GB ) ;

  /* Store result and return without error code */
  *result = Beta0 * m[0] * gpt_c / (-gpt_qe * stdsqrt(1.0 - Beta0 * Beta0)); ;
  return( 0 ) ;
}
