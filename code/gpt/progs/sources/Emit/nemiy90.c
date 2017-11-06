/* nemiy90.c - Calculate 90% normalized rms emittance for y y' subspace in [m-rad] */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "gdfa.h"

int nemiy90_func( double *result )
{
  int i, num, tmpnum ;
  double *nmacro, *y, *By, *Bz, *G  ;
  double *yc, *ypc ;
  double yy, ypyp, yyp ;
  double emiy, *emiyi, alphay, betay, gammay ;

  if( gdfmgetarr( "y",  &y,  &num    ) ||  num <=2     ||
      gdfmgetarr( "nmacro", &nmacro, &tmpnum ) || tmpnum!=num ||
      gdfmgetarr( "By", &By, &tmpnum ) ||  tmpnum!=num ||
      gdfmgetarr( "Bz", &Bz, &tmpnum ) ||  tmpnum!=num ||
      gdfmgetarr( "G" , &G,  &tmpnum ) ||  tmpnum!=num ) return(1) ;

  /* Allocate memory for new items xc, xpc */
  yc   = (double *)malloc( num*sizeof(double)) ;
  ypc  = (double *)malloc( num*sizeof(double)) ;
  emiyi= (double *)malloc( num*sizeof(double)) ;
  if( yc==NULL || ypc==NULL || emiyi==NULL ) { fprintf( stderr, "Not enough memory\n") ; return (1) ; }

  gdfasubavg(nmacro,yc, y, num) ;
  gdfasubavg(nmacro,ypc,By,num) ;

  yy   = gdfamean2(nmacro,yc ,yc ,num) ;
  ypyp = gdfamean2(nmacro,ypc,ypc,num) ;
  yyp  = gdfamean2(nmacro,yc ,ypc,num) ;

  emiy = stdsqrt( yy*ypyp-yyp*yyp) ;
  if( emiy==0.0 ) { *result=0.0 ; return(0) ; }

  alphay = -yyp/emiy ;
  betay  = yy/emiy ;
  gammay = ypyp/emiy ;

  for(i=0 ; i<num ; i++)
    emiyi[i] = gammay*yc[i]*yc[i] + 2*alphay*yc[i]*ypc[i] + betay*ypc[i]*ypc[i] ;
  qsort(emiyi,num,sizeof(*emiyi),gdfacmpdouble);

  *result = gdfamean(nmacro,G,num) * emiyi[(int)(0.90*num-0.5)] ;

  free( yc ) ;
  free( ypc ) ;
  free( emiyi ) ;

  return(0) ;
}
