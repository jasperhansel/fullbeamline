/* TMcyl010cavity.c - Cylindrical TM010 mode resonant cavity */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

#define x01 2.4048255576957727686216318

struct tmcav_info
{
  double d, radius ;
  double fac, phi ;
  double w010 ;
} ;

static int tmcav_sim(gptpar *par, double t, struct tmcav_info *info) ;


void TM010cylcavity_init(gptinit *init)
{
  struct tmcav_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=4 )
    gpterror( "Syntax: %s(ECS,d,R,const,phi)\n", gptgetname(init) ) ;

  info = (struct tmcav_info *)gptmalloc( sizeof(struct tmcav_info) ) ;

  info->d      = gptgetargdouble(init,1)/2 ;
  info->radius = gptgetargdouble(init,2) ;
  info->fac    = gptgetargdouble(init,3) ;
  info->phi    = gptgetargdouble(init,4) ;

  if( info->d<=0 || info->radius<=0 )
    gpterror( "Cavity dimensions must be larger than zero\n" ) ;

  if( info->fac==0 )
    gptwarning( "Const is equal to zero, no power in cavity\n" ) ;

  info->w010 = x01*gpt_c/info->radius ;

  gptaddEBelement( init, tmcav_sim, gptfree, GPTELEM_LOCAL, info ) ;
}


static int tmcav_sim(gptpar *par, double t, struct tmcav_info *info)
{
  double r ;
  double fac, radius, w010, phi ;
  double Bphi ;  

  if( fabs(Z) > info->d ) return( 0 ) ;
  r = sqrt(X*X + Y*Y) ;
  radius = info->radius ;
  if( r > radius ) return( 0 ) ;

  fac    = info->fac ;
  w010   = info->w010 ;
  phi    = info->phi ;

  EZ   = fac*bessj0(x01*r/radius)*sin(w010*t + phi) ;
  Bphi = fac*bessj1(x01*r/radius)/gpt_c*cos(w010*t + phi) ;

  gptrphi2carth(0,Bphi,X,Y,&BX,&BY) ;

  return( 1 ) ;
}
