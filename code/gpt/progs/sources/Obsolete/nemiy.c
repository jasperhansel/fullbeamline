/* nemiy.c - Calculate normalized emittance for x vx subspace */

/* written by Kees van der Geer */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "gdfa.h"

#define c 299792458.0

static void tocenter(   double *dest, double *source, int numx ) ;
static double rotation( double *x, double *y, double *vx, double *vy, int numx ) ;
static void unrotate(   double *x, double *y, double *vx, double *vy,
                        double omega, int numx ) ;
static double mean(     double *G, int numx ) ;
static double emit(     double *x, double *vx, int numx) ;
static double meansq(   double *a, double *b,  int numx) ;


int nemiy_func( double *result )
{
  int i, numx, numBx, numy, numBy, numG ;
  double *x,  *y,  *Bx,  *By, *G  ;
  double *xc, *yc, *vxc, *vyc ;
  double Gamma, omega, epsx, epsy   ;

  if( gdfmgetarr( "x" , &x,  &numx  ) ||  numx <=2    ||
      gdfmgetarr( "Bx", &Bx, &numBx ) ||  numBx!=numx ||
      gdfmgetarr( "y" , &y,  &numy  ) ||  numy !=numx ||
      gdfmgetarr( "By", &By, &numBy ) ||  numBy!=numx ||
      gdfmgetarr( "G" , &G,  &numG  ) ||  numG !=numx    ) return(1) ;

  /* Allocate memory for new items xc, yc, vxc, vyc
   * These are the manipulated values.
   */

   xc = (double *)malloc( numx*sizeof(double)) ;
   yc = (double *)malloc( numx*sizeof(double)) ;
   vxc= (double *)malloc( numx*sizeof(double)) ;
   vyc= (double *)malloc( numx*sizeof(double)) ;

   if (xc==NULL || yc==NULL || vxc==NULL || vyc==NULL)  
   {
     printf( "Not enough memory\n") ; return (1) ;
   }

   tocenter( xc,  x,  numx ) ;
   tocenter( yc,  y,  numx ) ;
   tocenter( vxc, Bx, numx ) ; 
   tocenter( vyc, By, numx ) ;
   
   /* change beta to velocity */
   for( i=0; i<numx; i++) { vxc[i] *= c ;   vyc[i] *= c ; } ;
       
   omega = rotation(xc, yc, vxc, vyc, numx) ;
   unrotate(        xc, yc, vxc, vyc, omega, numx) ;
   
   Gamma = mean(G, numx) ;
   epsx  = emit(xc, vxc, numx ) * Gamma/c ; /* stdsqrt(Gamma*Gamma-1) ; */
   epsy  = emit(yc, vyc, numx ) * Gamma/c ; /* stdsqrt(Gamma*Gamma-1) ; */
   *result = epsy ;
   /* let op, we hebben met gamma/c vermigvuldigd omdat we vergeten waren
      te werken met de hoek i.p.v. met de snelheid. Dus alle vx moet eigenlijk
      vx/vz worden
    */

   free( xc ); free( vxc ); free( yc ); free( vyc ) ;

   return(0) ;
}


static void tocenter( double *dest, double *source, int numx )
{
   int i ;
   double tmp=0.0 ;

   for(i=0 ; i<numx ; i++) tmp+=source[i] ;
   for(i=0 ; i<numx ; i++) dest[i] = source[i] - tmp/numx ;
   return ;
}

static double rotation(double *x, double *y, double *vx, double *vy, int numx )
{
  int i ;
  double tel, noem, r, vphi, omega ;

  tel = noem = 0 ;
  for(i=0 ; i<numx ; i++)
  {
    r     = stdsqrt( x[i]*x[i] + y[i]*y[i]) + 1.0e-9 ;
    vphi  = vy[i]*x[i]/r - vx[i]*y[i]/r ;
    tel  += vphi*r ;
    noem += r*r ;
  } ;
    omega = tel/noem ;
  return omega ;
}

static void unrotate( double *x,  double *y,  double *vx, double *vy,  
                      double omega, int numx )
{
   int i ;

   for(i=0 ; i<numx ; i++) vx[i] += omega*y[i] ;
   for(i=0 ; i<numx ; i++) vy[i] -= omega*x[i] ;
   return ;
}


static double mean( double *G, int numx )
{
   int i ;
   double tmp=0.0 ;

   for(i=0 ; i<numx ; i++) tmp+=G[i] ;
   return tmp/numx ;
}

static double emit(double *x, double *vx, int numx)
{
   double xx, yy, xy ;

   xx = meansq( x, x,numx) ;
   yy = meansq(vx,vx,numx) ;
   xy = meansq( x,vx,numx) ;

   return (4*stdsqrt(xx*yy-xy*xy)) ;
}

static double meansq( double *a, double *b, int numx )
{
   int i ;
   double temp = 0.0 ;

   for(i=0 ; i<numx ; i++) temp += a[i]*b[i] ;
   return temp/numx ;
}

