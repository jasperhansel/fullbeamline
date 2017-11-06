/* setbyemi.c - Set transverse y-emittance */

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

void setGByemittance_init(gptinit *init)
{
  gptparset *set ;
  gptinitpar *par ;
  char *name ;
  int i, num, mags ;
  double nemiy, oldnemiy, scale ;
  double *yc, *ypc, yy, ypyp, yyp ;

  if( gptgetargnum(init)<2 )
    gpterror( "Syntax: %s(set,emittance)\n", gptgetname(init) ) ;

  name  = gptgetargstring(init,1) ;
  nemiy = gptgetargdouble(init,2) ;

  /* Get particle set */
  if( gpttestparset( name )==NULL )
    gptwarning( "The particle set \"%s\" does not exist\n", name ) ;
  set = gptgetparset( name ) ;
  par = gptgetparsetpars( set,&num ) ;
  if( num<3 ) gpterror( "Insufficient number of particles to scale emittance\n" ) ;

  /* Calculate emittance */
  yc  = (double *)gptmalloc( num*sizeof(double)) ;
  ypc = (double *)gptmalloc( num*sizeof(double)) ;

  for(i=0 ; i<num ; i++) ypc[i] = par[i].GBr[1] ; /* yprime * Bz * G */
  for(i=0 ; i<num ; i++) yc[i]  = par[i].Wr[1]  ; /* y-coordinate */
  tocenter(yc, yc, num) ;
  tocenter(ypc,ypc,num) ;

  yy   = meansq(yc ,yc ,num) ;
  ypyp = meansq(ypc,ypc,num) ;
  yyp  = meansq(yc ,ypc,num) ;
  oldnemiy = stdsqrt( yy*ypyp-yyp*yyp) ;

  gptfree( yc ) ;
  gptfree( ypc ) ;

  /* Test scaling factor */
  if( oldnemiy==0.0 ) gpterror( "Initial imittance is zero. Scaling impossible\n" ) ;
  scale = nemiy / oldnemiy ;
  mags  = (int)fabs(log10(scale)) ;
  if( mags!=0 )
    gptwarning( "Scaling by over %d order%s of magnitide\n", mags, mags!=1 ? "s" : "" ) ;

  /* Scale particle distribution */
  for( i=0 ; i<num ; i++ ) par[i].GBr[1] *= scale ;
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
