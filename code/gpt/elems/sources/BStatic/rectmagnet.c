/* rectmagnet.c - Rectangular magnet with fringe fields */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 *
 * With special thanks to Bruno Muratory from Daresbury Laboratory
 * for providing us with the fringe field expressions
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

struct rectmagnet_info
{
  double a,b,Bo ;
  double dl,b1,b2 ;
} ;

static int rectmagnet_sim(gptpar *par,double t,struct rectmagnet_info *info) ;
static int rectmagnet_simple_sim(gptpar *par,double t,struct rectmagnet_info *info) ;

#define NGAP 10 /* Transient region is NGAP/b1 */

void rectmagnet_init(gptinit *init)
{
  struct rectmagnet_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=6 )
    gpterror( "Syntax: %s(ECS,a,b,Bfield,dl,b1,b2)\n", gptgetname(init) ) ;

  info = (struct rectmagnet_info *)gptmalloc( sizeof(struct rectmagnet_info) ) ;

  info->a  = gptgetargdouble(init,1)/2 ;
  info->b  = gptgetargdouble(init,2)/2 ;
  info->Bo = gptgetargdouble(init,3) ;
  info->dl = gptgetargdouble(init,4) ;
  info->b1 = gptgetargdouble(init,5) ;
  info->b2 = gptgetargdouble(init,6) ;

  if( info->b1<0 )
    gpterror( "%s: b1 must be nonnegative\n", gptgetname(init) ) ;

  if( info->b1==0 )
    gptaddEBelement( init, rectmagnet_simple_sim, gptfree, GPTELEM_LOCAL, info ) ;
  else
    gptaddEBelement( init, rectmagnet_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}


static int rectmagnet_sim(gptpar *par,double t,struct rectmagnet_info *info)
{
  double d ;    	/* Distance from edge */
  double f, h, Bn ;

  if( fabs(X)>info->a ) return( 0 ) ;

  d = fabs(Z)-info->b-info->dl ;

  if( d*info->b1 > NGAP ) return( 0 ) ;

  if( fabs(Y)*info->b1>=gpt_pi && d<=0) { gptremoveparticle(par) ; return(1) ; }

  f = info->b1*d + info->b2*(d*d-Y*Y) ;
  h = Y*(info->b1+2*info->b2*d) ;

  Bn = info->Bo/(1+2*exp(f)*cos(h)+exp(2*f)) ;
  BY = Bn*(1+exp(f)*cos(h)) ;
  BZ = -Bn*exp(f)*sin(h) ;

  if(Z<0) BZ=-BZ ;

  return(1) ;
}


static int rectmagnet_simple_sim(gptpar *par,double t,struct rectmagnet_info *info)
{
  if( fabs(X)>info->a ||
      fabs(Z)>info->b ) return( 0 ) ;

  BY = info->Bo ;

  return(1) ;
}
