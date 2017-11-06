/* gminmax.c: Remove particle outside energy interval */

#include <stdio.h>
#include <math.h>
#include "elem.h"

struct Gminmax_info
{
  double L ;
  double Gmin ;
  double Gmax ;
} ;

static int Gminmax_sim(gptpar *par,double t,struct Gminmax_info *info) ;

void Gminmax_init(gptinit *init)
{
  struct Gminmax_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=3 )
    gpterror( "Syntax: %s(ECS,L,Gmin,Gmax)\n", gptgetname(init) ) ;

  info = (struct Gminmax_info *)gptmalloc( sizeof(struct Gminmax_info) ) ;

  info->L        = gptgetargdouble(init,1) ;
  info->Gmin     = gptgetargdouble(init,2) ;
  info->Gmax     = gptgetargdouble(init,3) ;

  gptaddEBelement( init, Gminmax_sim, gptfree, GPTELEM_LOCAL, info ) ;
}


static int Gminmax_sim(gptpar *par,double t,struct Gminmax_info *info)
{
  if( fabs(Z)>info->L/2 ) return( 0 ) ;

  if( par->G < info->Gmin || par->G > info->Gmax ) gptremoveparticle(par) ;

  return( 1 ) ;
}
