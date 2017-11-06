/* nemirrms.c - Calculate normalized rms emittance for x x' y y' subspace in [m-rad] */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "gdfa.h"

int nemirrms_func( double *result )
{
  int num, tmpnum ;
  double *nmacro, *x,  *Bx,  *y, *By, *Bz, *G  ;
  double *xc, *xpc, *yc, *ypc ;
  double xx, xpxp, xxp, yy, ypyp, yyp, xy, xpyp, xpy, xyp ;

  if( gdfmgetarr( "x" , &x,  &num    ) ||  num <=2     ||
      gdfmgetarr( "nmacro", &nmacro, &tmpnum ) || tmpnum!=num ||
      gdfmgetarr( "Bx", &Bx, &tmpnum ) ||  tmpnum!=num ||
      gdfmgetarr( "y",  &y,  &tmpnum ) ||  tmpnum!=num ||
      gdfmgetarr( "By", &By, &tmpnum ) ||  tmpnum!=num ||
      gdfmgetarr( "Bz", &Bz, &tmpnum ) ||  tmpnum!=num ||
      gdfmgetarr( "G" , &G,  &tmpnum ) ||  tmpnum!=num ) return(1) ;

  /* Allocate memory for new items xc, xpc */
  xc = (double *)malloc( num*sizeof(double)) ;
  xpc= (double *)malloc( num*sizeof(double)) ;
  yc = (double *)malloc( num*sizeof(double)) ;
  ypc= (double *)malloc( num*sizeof(double)) ;
  if (xc==NULL || xpc==NULL || yc==NULL || ypc==NULL )  { fprintf( stderr, "Not enough memory\n") ; return (1) ; }

  gdfasubavg(nmacro, xc, x, num ) ;
  gdfasubavg(nmacro, xpc,Bx,num ) ;
  gdfasubavg(nmacro, yc, y, num ) ;
  gdfasubavg(nmacro, ypc,By,num ) ;

  xx   = gdfamean2(nmacro,xc ,xc ,num) ;
  xpxp = gdfamean2(nmacro,xpc,xpc,num) ;
  xxp  = gdfamean2(nmacro,xc ,xpc,num) ;

  yy   = gdfamean2(nmacro,yc ,yc ,num) ;
  ypyp = gdfamean2(nmacro,ypc,ypc,num) ;
  yyp  = gdfamean2(nmacro,yc ,ypc,num) ;

  xy   = gdfamean2(nmacro,xc ,yc ,num) ;
  xpyp = gdfamean2(nmacro,xpc,ypc,num) ;
  xyp  = gdfamean2(nmacro,xc ,ypc,num) ;
  xpy  = gdfamean2(nmacro,xpc,yc ,num) ;

  *result = stdsqrt( (xx*xpxp-xxp*xxp)*(yy*ypyp-yyp*yyp) ) - fabs(xy*xpyp-xyp*xpy) ;
  *result = gdfamean(nmacro,G,num) * stdsqrt( *result ) ;

  free( xc ) ;
  free( xpc ) ;
  free( yc ) ;
  free( ypc ) ;

  return(0) ;
}
