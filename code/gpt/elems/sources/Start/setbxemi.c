/* setbxemi.c - Set transverse x-emittance */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "elem.h"

/* Helper functions */
static void tocenter( double *dest, double *source, int num ) ;
static double meansq( double *a, double *b, int num ) ;
static double stdsqrt(double x) ;

void setGBxemittance_init(gptinit *init)
{
  gptparset *set ;
  gptinitpar *par ;
  char *name ;
  int i, num, mags ;
  double nemix, oldnemix, scale ;
  double *xc, *xpc, xx, xpxp, xxp ;

  if( gptgetargnum(init)<2 )
    gpterror( "Syntax: %s(set,emittance)\n", gptgetname(init) ) ;

  name  = gptgetargstring(init,1) ;
  nemix = gptgetargdouble(init,2) ;

  /* Get particle set */
  if( gpttestparset( name )==NULL )
    gptwarning( "The particle set \"%s\" does not exist\n", name ) ;
  set = gptgetparset( name ) ;
  par = gptgetparsetpars( set,&num ) ;
  if( num<3 ) gpterror( "Insufficient number of particles to scale emittance\n" ) ;

  /* Calculate emittance */
  xc  = (double *)gptmalloc( num*sizeof(double)) ;
  xpc = (double *)gptmalloc( num*sizeof(double)) ;

  for(i=0 ; i<num ; i++) xpc[i] = par[i].GBr[0] ; /* xprime * Bz * G */
  for(i=0 ; i<num ; i++) xc[i]  = par[i].Wr[0]  ; /* x-coordinate */
  tocenter(xc, xc, num) ;
  tocenter(xpc,xpc,num) ;

  xx   = meansq(xc ,xc ,num) ;
  xpxp = meansq(xpc,xpc,num) ;
  xxp  = meansq(xc ,xpc,num) ;
  oldnemix = stdsqrt( xx*xpxp-xxp*xxp) ;

  gptfree( xc ) ;
  gptfree( xpc ) ;

  /* Test scaling factor */
  if( oldnemix==0.0 ) gpterror( "Initial imittance is zero. Scaling impossible\n" ) ;
  scale = nemix / oldnemix ;
  mags  = (int)fabs(log10(scale)) ;
  if( mags!=0 )
    gptwarning( "Scaling by over %d order%s of magnitide\n", mags, mags!=1 ? "s" : "" ) ;

  /* Scale particle distribution */
  for( i=0 ; i<num ; i++ ) par[i].GBr[0] *= scale ;
}

static void tocenter( double *dest, double *source, int num )
{
  int i ;
  double tmp=0.0 ;

  for(i=0 ; i<num ; i++) tmp+=source[i] ; tmp /= num ;
  for(i=0 ; i<num ; i++) dest[i] = source[i] - tmp;
  return ;
}


static double meansq( double *a, double *b, int num )
{
  int i ;
  double tmp = 0.0 ;

  for(i=0 ; i<num ; i++) tmp += a[i]*b[i] ;
  return tmp/num ;
}


static double stdsqrt(double x)
{
  if( x<0.0 ) return( 0.0 ) ;
  return( sqrt(x) ) ;
}
