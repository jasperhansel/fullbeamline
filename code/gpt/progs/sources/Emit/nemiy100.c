/* nemiy100.c - Calculate 100% normalized rms emittance for y y' subspace in [m-rad] */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "gdfa.h"

int nemiy100_func( double *result )
{
  int i, num, tmpnum ;
  double *nmacro, *y, *By, *Bz, *G  ;
  double *yc, *ypc ;
  double yy, ypyp, yyp ;
  double emiy, emiy100, emiyi, alphay, betay, gammay ;

  if( gdfmgetarr( "y",  &y,  &num    ) ||  num <=2     ||
      gdfmgetarr( "nmacro", &nmacro, &tmpnum ) || tmpnum!=num ||
      gdfmgetarr( "By", &By, &tmpnum ) ||  tmpnum!=num ||
      gdfmgetarr( "Bz", &Bz, &tmpnum ) ||  tmpnum!=num ||
      gdfmgetarr( "G" , &G,  &tmpnum ) ||  tmpnum!=num ) return(1) ;

  /* Allocate memory for new items xc, xpc */
  yc = (double *)malloc( num*sizeof(double)) ;
  ypc= (double *)malloc( num*sizeof(double)) ;
  if( yc==NULL || ypc==NULL ) { fprintf( stderr, "Not enough memory\n") ; return (1) ; }

  gdfasubavg(nmacro,yc, y,  num) ;
  gdfasubavg(nmacro,ypc,By,num) ;

  yy   = gdfamean2(nmacro,yc ,yc ,num) ;
  ypyp = gdfamean2(nmacro,ypc,ypc,num) ;
  yyp  = gdfamean2(nmacro,yc ,ypc,num) ;

  emiy = stdsqrt( yy*ypyp-yyp*yyp) ;
  if( emiy==0.0 ) { *result=0.0 ; return(0) ; }

  alphay = -yyp/emiy ;
  betay  = yy/emiy ;
  gammay = ypyp/emiy ;

  emiy100 = 0.0 ;
  for(i=0 ; i<num ; i++)
  {
    emiyi = gammay*yc[i]*yc[i] + 2*alphay*yc[i]*ypc[i] + betay*ypc[i]*ypc[i] ;
    if( emiyi>emiy100 ) emiy100 = emiyi ;
  }
  *result = gdfamean(nmacro,G,num) * emiy100 ;

  free( yc ) ;
  free( ypc ) ;

  return(0) ;
}
