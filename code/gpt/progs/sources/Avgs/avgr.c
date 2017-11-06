/* avgr.c - calculate average r */

/* written by Kees van der Geer */

#include <math.h>

#include "gdfa.h"

int avgr_func( double *result )
{
  /* Declarations */
  int i, num, tmpnum ;
  double avgx, avgy ;
  double r, sumr=0, sumnr=0, sumn=0 ;
  double *x, *y, *nmacro ;

  /* Get selected arrays from GDFA kernel */
  if( gdfmgetarr( "x", &x, &num ) || num<1  ||
      gdfmgetarr( "y", &y, &tmpnum ) || tmpnum!=num ||
      gdfmgetarr( "nmacro", &nmacro, &tmpnum ) || tmpnum!=num ) return(1) ;

  avgx = gdfamean(nmacro,x,num) ;
  avgy = gdfamean(nmacro,y,num) ;
  for(i=0 ; i<num ; i++)
  {
    r = sqrt((x[i]-avgx)*(x[i]-avgx)+(y[i]-avgy)*(y[i]-avgy)) ;
    sumr  += r ;
    sumnr += nmacro[i]*r ;
    sumn  += nmacro[i] ;
  }

  /* Store result and return without error code */
  if( sumn==0 )
    *result = sumr/num ;
  else
    *result = sumnr/sumn ;

  return( 0 ) ;
}
