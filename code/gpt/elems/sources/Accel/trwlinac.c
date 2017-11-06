/* trwlinac.c - travelling wave linac without beamloading */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 *
 * written by Marieke de Loos after twlinac.c by Kees van der Geer
 */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "elem.h"

#define SHAPELEN 0.001        /* Shaping over 1 mm */

struct trwlinac_info
{
  double ao ;	  /* initial attenuation constant			*/
  double Go ;	  /* design gamma at entrance                           */
  double phi ;	  /* rf phase						*/
  double w  ;     /* radial frequency					*/
  double L ;      /* half length				        */

  double E1 ;     /* constant gradient field component			*/
  double Fomc2 ;  /* qE/(mcc)					        */
  double BGo ;    /* beta*gamma for design values                       */
  double ko ;     /* wavenumber                                         */
} ;

static int trwlinac_sim(gptpar *par,double t,struct trwlinac_info *info) ;


void trwlinac_init(gptinit *init)
{
  struct trwlinac_info *info ;
  double ao, rs, Po, P, Go, w, Eo ;
  double thetao ; /* bunch phase w.r.t. wave, design value              */

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=9 )
    gpterror( "Syntax: %s(ECS,ao,rs,Po,P,Go,thetao,phi,w,L)\n", gptgetname(init) ) ;
  
  info = (struct trwlinac_info *)gptmalloc( sizeof(struct trwlinac_info) ) ;

  ao    = info->ao    = gptgetargdouble(init,1) ;
  rs    =               gptgetargdouble(init,2) ;
  Po    =               gptgetargdouble(init,3) ;
  P     =               gptgetargdouble(init,4) ;
  Go    = info->Go    = gptgetargdouble(init,5) ;
  thetao              = gptgetargdouble(init,6) ;
          info->phi   = gptgetargdouble(init,7) ;
  w     = info->w     = gptgetargdouble(init,8) ;
          info->L     = gptgetargdouble(init,9)/2 ;

  Eo          = sqrt(2*ao*rs*Po) * cos(thetao) ;     /*design value */ 
  info->Fomc2 = -Eo * gpt_qe/(gpt_me*gpt_c*gpt_c) ;  /* idem        */
  info->BGo   = sqrt(Go*Go - 1) ; 
  info->E1    = sqrt(2*ao*rs*P) ;                    /*real top value*/
  info->ko    = w/gpt_c ;

  gptaddEBelement( init, trwlinac_sim, gptfree, GPTELEM_LOCAL, info ) ;
}


static int trwlinac_sim(gptpar *par,double t,struct trwlinac_info *info)
{
  double ao, Go, phi, w, L, Fomc2, BGo ;
  double theta          ;     /* actual phase                          */
  double E1             ;     /* actual Eo                             */
  double Er, Bp	 	;     /* field to be returned		       */
  double intkz	 	;     /* integral kz from 0 to z	       */
  double beta2 	 	;     /* beta(z)^2 for unloaded linac	       */
  double Zp, r	  	;     /* particle coordinates		       */
  double ko, kt, kz	;     /* wavenumbers                           */
  double shape 	        ;     /* smoothing factor for border	       */

  ao	= info->ao    ;
  Go    = info->Go    ;
  phi   = info->phi   ;
  w     = info->w     ;
  L     = info->L     ;       

  E1    = info->E1    ;
  Fomc2 = info->Fomc2 ;
  BGo   = info->BGo   ;
  ko    = info->ko    ;

  Zp    = Z + L       ;       /* Position relative to begin of linac  */
   
  /* Shaping */
  if( fabs(Z) > L+3*SHAPELEN ) return( 0 ) ; 
  shape = erfshape((Z+L)/SHAPELEN) * erfshape((L-Z)/SHAPELEN) ;

  if ((1.0-2.0*ao*Zp) <= 0.0 ) gpterror( "structure too long\n");

  intkz  = (ko/Fomc2) * (sqrt( SQR(Fomc2*Zp + Go) - 1.0) - BGo) ;
  theta  = w*t + phi - intkz ;
  
  beta2  = 1.0 - 1.0/SQR(Go + Fomc2*Zp ) ;
  if( beta2<=0.0 ) gpterror( "beta^2 is negative or zero!\n" ) ;
  kz     = ko/sqrt(beta2) ; 
  kt     = SQR(kz)-SQR(ko) ;
  if( kt<=0.0 ) gpterror( "kz^2-ko^2 is negative or zero!\n" ) ;
  kt     = sqrt(kt) ;

  r  = sqrt( SQR(X)+SQR(Y) ) ;   

  EZ = shape * E1*sin(theta)*bessi0(kt*r) ;
  Er = shape * E1*cos(theta)*bessi1(kt*r) * kz/kt ;
  Bp = Er/gpt_c ;

  /* Convert to carthesian coordinate system. Could be improved. */
  if( r>DBL_EPSILON )
  {
    EX = X * Er / r ;
    EY = Y * Er / r ;
    BX =-Y * Bp / r ;
    BY = X * Bp / r ;
  }

  return( 1 ) ;
}
