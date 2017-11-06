/* ezcell.c - Tabulated Ez-file, using spline */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

struct ezcell_info
{
  unsigned int n ;
  double *z ;
  double *ezcell ;
  double *ezcell2 ;
  double phi ;
  double omega ;
} ;

static int ezcell_sim(gptpar *par,double t,struct ezcell_info *info) ;
static void ezcell_exit( struct ezcell_info *info ) ;


void ezcell_init(gptinit *init)
{
  struct ezcell_info *info ;
  double Ezef ;
  unsigned int i ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=4 )
    gpterror( "Syntax: %s(ECS,filename,Ezef,phi,w)\n", gptgetname(init) ) ;

  info = (struct ezcell_info *)gptmalloc( sizeof(struct ezcell_info) ) ;

  Ezef=gptgetargdouble(init,2) ;

  readtabfile( init, gptgetargstring(init,1), &info->z, &info->ezcell, &info->n ) ;

  info->ezcell2 = (double *)gptmalloc(info->n*sizeof(double)) ;

  for(i=0 ; i<info->n ; i++) info->ezcell[i] *= Ezef ;

  if( spline( info->z, info->ezcell, info->ezcell2, info->n, NATSPLINE, NATSPLINE ) )
    gpterror( "Error calculating spline\n" ) ;

  info->phi = gptgetargdouble(init,3) ;
  info->omega=gptgetargdouble(init,4) ;

  gptaddEBelement( init, ezcell_sim, ezcell_exit, GPTELEM_LOCAL, info ) ;
}


static int ezcell_sim(gptpar *par,double t,struct ezcell_info *info)
{
  double y, dy ;
  double omega, phase, cosphase, sinphase ;
  double Ez, Eror, Bphior ;

  if( splintd( info->z, info->ezcell, info->ezcell2, info->n, Z, &y, &dy ) ) return( 0 ) ;

  omega = info->omega ;
  phase = omega*t+info->phi ;
  cosphase = cos(phase) ;
  sinphase = sin(phase) ;

  Ez     =  cosphase*y ;
  Eror   = -cosphase*dy*0.5 ;
  Bphior = -sinphase*y*omega/(2.0*gpt_c*gpt_c) ;

  EX =  X*Eror ;
  EY =  Y*Eror ;
  EZ =  Ez ;
  BX = -Y*Bphior ;
  BY =  X*Bphior ;

  return( 1 ) ;
}


static void ezcell_exit( struct ezcell_info *info )
{
  gptfree( info->z ) ;
  gptfree( info->ezcell ) ;
  gptfree( info->ezcell2 ) ;
  gptfree( info ) ;
}
