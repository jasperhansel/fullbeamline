/* trwlinbm.c - travelling wave linac with beamloading modelled by ODE's */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 *
 * written by Marieke de Loos
 */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "elem.h"

#define SHAPELEN 0.001        /* Shaping over 1 mm */

struct trwlinbm_info
{
  double ao ;	  /* initial attenuation constant			*/
  double rs ;	  /* shunt impedance					*/
  double Ib ;	  /* beam current					*/
  double Go ;	  /* design gamma at entrance				*/
  double phi;	  /* rf phase						*/
  double w  ;     /* radial frequency					*/
  double L ;      /* half length				        */

  double E1 ;     /* constant gradient field component			*/
  double Fomc2 ;  /* qE/(mcc)					        */
  double BGo ;    /* beta*gamma for design values                       */
  double ko ;     /* wavenumber                                         */

  int offset ;    /* Position of u and v in ODElist */
  double *ODEx ;  /* Filled in by trwlinbm_fpr */
  double *ODEp ;  /* Filled in by trwlinbm_fpr */
} ;

static void trwlinbm_fpr(double t, double *x, double *p, void *vinfo) ;
static int  trwlinbm_sim(gptpar *par,double t,struct trwlinbm_info *info) ;
static void trwlinbm_out(double t, double *dt, double *x, struct trwlinbm_info *info ) ;


void trwlinbm_init(gptinit *init)
{
  struct trwlinbm_info *info ;
  double ao, rs, Po, P, Go, thetao, w, Eo, u, v ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=10 )
    gpterror( "Syntax: %s(ECS,ao,rs,Po,P,Ib,Go,thetao,phi,w,L)\n", gptgetname(init) ) ;
  
  info = (struct trwlinbm_info *)gptmalloc( sizeof(struct trwlinbm_info) ) ;

  ao = info->ao    = gptgetargdouble(init,1) ;
  rs = info->rs    = gptgetargdouble(init,2) ;
  Po =               gptgetargdouble(init,3) ;
  P  =               gptgetargdouble(init,4) ;
       info->Ib    = gptgetargdouble(init,5) ;
  Go = info->Go    = gptgetargdouble(init,6) ;
  thetao           = gptgetargdouble(init,7) ;
       info->phi   = gptgetargdouble(init,8) ;
  w  = info->w     = gptgetargdouble(init,9)  ;
       info->L     = gptgetargdouble(init,10)/2 ;

  Eo          = sqrt(2*ao*rs*Po) * cos(thetao) ;     /*design value */ 
  info->Fomc2 = -Eo * gpt_qe/(gpt_me*gpt_c*gpt_c) ;  /* idem        */
  info->BGo   = sqrt(Go*Go - 1) ; 
  info->E1    = sqrt(2*ao*rs*P) ;                    /*real top value*/
  info->ko    = w/gpt_c ;

  u = 0 ;
  v = 0 ;

/* Register function to kernel */
  gptaddEBelement( init, trwlinbm_sim, gptfree, GPTELEM_LOCAL, info ) ;
  toutaddfnc( trwlinbm_out, info ) ;
  odeaddvar( &odeinfo, 2, &info->offset ) ;
  odeinfo.x[info->offset+0] = u ;
  odeinfo.x[info->offset+1] = v ;
  odeaddfprfunction( ODEFNC_INI,trwlinbm_fpr,info ) ;
}


static void trwlinbm_fpr(double t, double *x, double *p, void *vinfo)
{
  struct trwlinbm_info *info = (struct trwlinbm_info *)vinfo ;

  info->ODEx = x ;
  info->ODEp = p ;

  p[info->offset+0] = 0 ;
  p[info->offset+1] = 0 ;
}


static int trwlinbm_sim(gptpar *par,double t,struct trwlinbm_info *info)
{
  double ao, rs, Ib, Go, phi, w, L, Fomc2, BGo ;
  double theta          ;     /* actual phase                          */
  double E1             ;     /* actual Eo                             */
  double Er, Bp	 	;     /* field to be returned		       */
  double intkz	 	;     /* integral kz from 0 to z	       */
  double beta2 	 	;     /* beta(z)^2 for unloaded linac	       */
  double Zp, r	  	;     /* particle coordinates		       */
  double ko, kt, kz	;     /* wavenumbers                           */
  double shape 	        ;     /* smoothing factor for border	       */
  double u, v, fac      ;     /* ODE parameters                        */
  double sintheta, costheta, vz ;

  ao	= info->ao    ;
  rs	= info->rs    ;
  Ib	= info->Ib    ;
  Go    = info->Go    ;
  phi   = info->phi   ;
  w     = info->w     ;
  L     = info->L     ;              

  E1    = info->E1    ;
  Fomc2 = info->Fomc2 ;
  BGo   = info->BGo   ;
  ko    = info->ko    ;

  u     = info->ODEx[info->offset+0] ;
  v     = info->ODEx[info->offset+1] ;

  Zp    = Z + L       ;              /* Relative to begin of linac */
   
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

  sintheta = sin(theta) ;
  costheta = cos(theta) ;

  EZ = shape * (E1*sintheta + u*sintheta + v*costheta)*bessi0(kt*r) ;
  Er = shape * (E1*costheta + u*costheta - v*sintheta)*bessi1(kt*r) * kz/kt ;
  Bp = Er/gpt_c ;

  /* Convert to carthesian coordinate system. Could be improved. */
  if( r>DBL_EPSILON )
  {
    EX = X * Er / r ;
    EY = Y * Er / r ;
    BX =-Y * Bp / r ;
    BY = X * Bp / r ;
  }

  vz   = gpt_c*par->GBr[2]/par->G ;
  fac  = (ao*Ib*rs)/(1-2*ao*Zp)*vz/numpar ;

  info->ODEp[info->offset+0] += fac*sintheta ;
  info->ODEp[info->offset+1] += fac*costheta ;

  return( 1 ) ;
}


static void trwlinbm_out ( double t, double *dt, double *x, struct trwlinbm_info *info )
{
  double u, v, Eb, phib  ;

  u  = x[info->offset+0] ;
  v  = x[info->offset+1] ;

  Eb   = sqrt(u*u+v*v) ; 
  phib = (Eb>0) ? atan2(v,u) : 0 ;

  gpswval( "Eb", Eb ) ;
  gpswval( "phib", phib ) ;
}
