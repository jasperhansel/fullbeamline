/* CSgammay.c - Calculate Courant-Snyder gamma for y y' subspace in [rad/m] */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "gdfa.h"

int CSgammay_func( double *result )
{
  int num, tmpnum ;
  double *nmacro, *y, *By, *Bz, *G  ;
  double *yc, *ypc ;
  double avgBz, yy, ypyp, yyp, emiy ;

  if( gdfmgetarr( "y",  &y,  &num    ) ||  num <=2     ||
      gdfmgetarr( "nmacro", &nmacro, &tmpnum ) || tmpnum!=num ||
      gdfmgetarr( "By", &By, &tmpnum ) ||  tmpnum!=num ||
      gdfmgetarr( "Bz", &Bz, &tmpnum ) ||  tmpnum!=num ||
      gdfmgetarr( "G" , &G,  &tmpnum ) ||  tmpnum!=num ) return(1) ;

  /* Allocate memory for new items xc, xpc */
  yc = (double *)malloc( num*sizeof(double)) ;
  ypc= (double *)malloc( num*sizeof(double)) ;
  if( yc==NULL || ypc==NULL ) { fprintf( stderr, "Not enough memory\n") ; return (1) ; }

  gdfasubavg(nmacro,yc, y, num) ;
  gdfasubavg(nmacro,ypc,By,num) ;

  avgBz= gdfamean (nmacro,Bz ,num) ;
  yy   = gdfamean2(nmacro,yc ,yc ,num) ;
  ypyp = gdfamean2(nmacro,ypc,ypc,num) ;
  yyp  = gdfamean2(nmacro,yc ,ypc,num) ;

  emiy = stdsqrt( yy*ypyp-yyp*yyp) ;
  if( emiy==0.0 || avgBz==0.0 ) return( 1 ) ;

  *result = ypyp/(emiy*avgBz) ;

  free( yc ) ;
  free( ypc ) ;

  return(0) ;
}
