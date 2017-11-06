/* nemiz90.c - Calculate 90% normalized rms emittance for t E subspace in [eV*s] */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "gdfa.h"

int nemiz90_func( double *result )
{
  int i, num, tmpnum ;
  double *nmacro, *z,  *Bz,  *G,  *t,  *m, m0;
  double *tc, *Gc ;
  double GG, tt, tG ;
  double emiz, *emizi, alphaz, betaz, gammaz ;

  /* Retreive G and m coordinates */
  if( gdfmgetarr( "G" , &G,  &num    ) ||  num<=2      ||
      gdfmgetarr( "nmacro", &nmacro, &tmpnum ) || tmpnum!=num ||
      gdfmgetarr( "m" , &m,  &tmpnum ) ||  tmpnum!=num ) return(1) ;

   /* check particles: */
   m0 = m[0];
   for( i=0 ; i<num ; i++) if( m[i]!=m0 )
     { fprintf( stderr, "Nemizrms works only for particles of the same type!\n") ; return(1) ; }

  /* Allocate memory for Gc and tc arrays */
  Gc   = (double *)malloc( num*sizeof(double)) ;
  tc   = (double *)malloc( num*sizeof(double)) ;
  emizi= (double *)malloc( num*sizeof(double)) ;

  if( Gc==NULL || tc==NULL || emizi==NULL ) 
    { fprintf( stderr, "Not enough memory\n") ; return (1) ; }

  /* Calculate Gc and tc arrays */
  gdfasubavg(nmacro,Gc,G,num) ;
  if( gdfmgetarr( "t", &t, &tmpnum ))
  {
    /* No time information present: Retrieve z and Bz coordinates */
    if( gdfmgetarr( "z" , &z,  &tmpnum ) ||  tmpnum!=num ||
	gdfmgetarr( "Bz", &Bz, &tmpnum ) ||  tmpnum!=num ) return(1) ;

    gdfasubavg(nmacro,tc,z,num) ; /* tc temporarily contains centered z coordinates! */
    for(i=0 ; i<num ; i++) if( Bz[i]!=0.0 )
      tc[i] /= - (Bz[i]*gpt_c) ;
    else 
      { fprintf( stderr, "Particle beta must not be zero\n" ) ; return( 1 ) ; }
    gdfasubavg(nmacro,tc,tc,num) ;
  } else	
  {
    /* Simple case: Time information */
    if( tmpnum!=num  ) return(1) ;
    gdfasubavg(nmacro,tc,t,num) ;
  }

  /* Calculate emittance */
  GG = gdfamean2(nmacro,Gc,Gc,num) ;
  tt = gdfamean2(nmacro,tc,tc,num) ;
  tG = gdfamean2(nmacro,tc,Gc,num) ;

  emiz = stdsqrt(GG*tt-tG*tG) ;
  if( emiz==0.0 ) { *result=0.0 ; return(0) ; }

  alphaz = -tG/emiz ;
  betaz  = tt/emiz ;
  gammaz = GG/emiz ;

  for(i=0 ; i<num ; i++)
    emizi[i] = gammaz*tc[i]*tc[i] + 2*alphaz*tc[i]*Gc[i] + betaz*Gc[i]*Gc[i] ;
  qsort(emizi,num,sizeof(*emizi),gdfacmpdouble);

  *result = emizi[(int)(0.90*num-0.5)] * m0*gpt_c*gpt_c / -gpt_qe;

  /* Free memory */
  free( tc ) ;
  free( Gc ) ;
  free( emizi ) ;

  return(0) ;
}
