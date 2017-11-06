/* bzsol.c - Bz solenoid */

/* BCM-General Particle Simulation: Bz-solenoid
 *
 * This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"


struct bzsolenoid_info
{
  double R2 ;
  double zlen ;
  double mu0nI ;
} ;

static int bzsolenoid_sim(gptpar *par,double t,struct bzsolenoid_info *info) ;


void bzsolenoid_init( gptinit *init )
{
  struct bzsolenoid_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=3 )
    gpterror( "Syntax: %s(ECS,R,L,nI)\n", gptgetname(init) ) ;

  info = (struct bzsolenoid_info *)gptmalloc( sizeof(struct bzsolenoid_info) ) ;

  info->R2 = SQR(gptgetargdouble(init,1)) ;
  info->zlen = gptgetargdouble(init,2)/2 ;
  info->mu0nI = gpt_mu0*gptgetargdouble(init,3) ;

  gptaddEBelement( init, bzsolenoid_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}


static int bzsolenoid_sim(gptpar *par,double t,struct bzsolenoid_info *info)
{
  double zpL, zmL, a, b, Br ;

  zpL = Z + info->zlen ;
  zmL = Z - info->zlen ;

  a = sqrt( zpL*zpL + info->R2 ) ;
  b = sqrt( zmL*zmL + info->R2 ) ;

  Br = -info->mu0nI*info->R2*( 1.0/(a*a*a) - 1.0/(b*b*b) ) * 0.25 ;

  BX = X * Br ; 
  BY = Y * Br ;
  BZ = info->mu0nI*( zpL/a - zmL/b ) * 0.5 ;

  return( 1 ) ;
}
