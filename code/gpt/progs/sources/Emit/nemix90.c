/* nemix90.c - Calculate 90% normalized rms emittance for x x' subspace in [m-rad] */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "gdfa.h"

int nemix90_func( double *result )
{
  int i, num, tmpnum ;
  double *nmacro, *x,  *Bx, *Bz, *G  ;
  double *xc, *xpc ;
  double xx, xpxp, xxp ;
  double emix, *emixi, alphax, betax, gammax ;

  if( gdfmgetarr( "x" , &x,  &num    ) ||  num <=2     ||
      gdfmgetarr( "nmacro", &nmacro, &tmpnum ) || tmpnum!=num ||
      gdfmgetarr( "Bx", &Bx, &tmpnum ) ||  tmpnum!=num ||
      gdfmgetarr( "Bz", &Bz, &tmpnum ) ||  tmpnum!=num ||
      gdfmgetarr( "G" , &G,  &tmpnum ) ||  tmpnum!=num ) return(1) ;

  /* Allocate memory for new items xc, xpc */
  xc   = (double *)malloc( num*sizeof(double)) ;
  xpc  = (double *)malloc( num*sizeof(double)) ;
  emixi= (double *)malloc( num*sizeof(double)) ;
  if( xc==NULL || xpc==NULL || emixi==NULL )  { fprintf( stderr, "Not enough memory\n") ; return (1) ; }

  gdfasubavg(nmacro,xc, x, num) ;
  gdfasubavg(nmacro,xpc,Bx,num) ;

  xx   = gdfamean2(nmacro,xc ,xc ,num) ;
  xpxp = gdfamean2(nmacro,xpc,xpc,num) ;
  xxp  = gdfamean2(nmacro,xc ,xpc,num) ;

  emix = stdsqrt( xx*xpxp-xxp*xxp) ;
  if( emix==0.0 ) { *result=0.0 ; return(0) ; }

  alphax = -xxp/emix ;
  betax  = xx/emix ;
  gammax = xpxp/emix ;

  for(i=0 ; i<num ; i++)
    emixi[i] = gammax*xc[i]*xc[i] + 2*alphax*xc[i]*xpc[i] + betax*xpc[i]*xpc[i] ;
  qsort(emixi,num,sizeof(*emixi),gdfacmpdouble);

  *result = gdfamean(nmacro,G,num) * emixi[(int)(0.90*num-0.5)] ;

  free( xc ) ;
  free( xpc ) ;
  free( emixi ) ;

  return(0) ;
}
