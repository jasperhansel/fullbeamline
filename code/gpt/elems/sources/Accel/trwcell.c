/* trwcell.c - Travelling wave resonant cavity */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "elem.h"

struct trwcell_info
{
  double M ;
  double phi ;
  double omega ;
  double zlen ;
  double Kz ;
  double Kt ;
} ;

static int trwcell_sim(gptpar *par,double t,struct trwcell_info *info) ;


void trwcell_init(gptinit *init)
{
  struct trwcell_info *info ;
  double tmp ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=4 )
    gpterror( "Syntax: %s(ECS,Ezef,phi,w,L)\n", gptgetname(init) ) ;

  info = (struct trwcell_info *)gptmalloc( sizeof(struct trwcell_info) ) ;

  info->M     = gptgetargdouble(init,1) ;
  info->phi   = gptgetargdouble(init,2) ;
  info->omega = gptgetargdouble(init,3) ;
  info->zlen  = gptgetargdouble(init,4)/2 ;

  info->M *= (3.0*sqrt(3.0))/(2.0*gpt_pi) ;
  info->Kz = gpt_pi/(3.0*info->zlen) ;

  tmp = SQR(info->Kz)-SQR(info->omega/gpt_c) ;
  if( tmp<=0 ) gpterror( "Kz^2-(omega/c)^2 is negative or zero!\n" ) ;
  info->Kt = sqrt(tmp) ;

  gptaddEBelement( init, trwcell_sim, gptfree, GPTELEM_LOCAL, info ) ;
}


static int trwcell_sim(gptpar *par,double t,struct trwcell_info *info)
{
  double phase ;
  double r ;
  double Ez ;
  double Er ;
  double Bp ;

  if( fabs(Z) > info->zlen ) return( 0 ) ;

  phase = info->omega*t + info->phi - info->Kz*Z ;
  r = sqrt( X*X+Y*Y ) ;

  Ez = info->M*sin(phase)*bessi0(info->Kt*r) ;
  Er = info->M*cos(phase)*bessi1(info->Kt*r) * info->Kz/info->Kt ;
  Bp = Er/gpt_c ;

  /* Convert to carthesian coordinate system. Could be improved. */
  if( r>DBL_EPSILON )
  {
    EX =  X * Er / r ;
    EY =  Y * Er / r ;
    BX = -Y * Bp / r ;
    BY =  X * Bp / r ;
  }

  EZ = Ez ;

  return( 1 ) ;
}
