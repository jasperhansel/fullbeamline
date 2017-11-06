/* multisli.c - Kill particle when it does not cross a slit opening */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "elem.h"

struct multislit_info
{
  double ap ;     /* a/(2 d)  */
  double b ;
  double L ;
  double d ;
  int N ;
} ;

static int multislit_sim(gptpar *par,double t,struct multislit_info *info) ;


void multislit_init(gptinit *init)
{
  struct multislit_info *info ;
  double a ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=5 )
    gpterror( "Syntax: %s(ECS,a,b,L,d,N)\n", gptgetname(init) ) ;

  info = (struct multislit_info *)gptmalloc( sizeof(struct multislit_info) ) ;

        a = gptgetargdouble(init,1)/2 ;
  info->b = gptgetargdouble(init,2)/2 ;
  info->L = gptgetargdouble(init,3)/2 ;
  info->d = gptgetargdouble(init,4) ;
  info->N = gptgetargint   (init,5) ;

  info->ap = a/info->d ;

  if( info->N<=0 ) gpterror("Number of slits must be positive\n") ; 
  if( 2*a>=info->d ) gpterror("Slit opening is larger than the distance between the slits\n") ; 

  gptaddEBelement( init, multislit_sim, gptfree, GPTELEM_LOCAL, info ) ;
}


static int multislit_sim(gptpar *par,double t,struct multislit_info *info)
{
  double d, xp, ap ;
  int N ;

  if( fabs(Z) > info->L ) return( 0 ) ;

  d = info->d ;
  N = info->N ;
  ap = info->ap ;
  xp = X/d - floor(X/d);                          /* Suggest: modf */

  if( fabs(Y)>info->b ||                          /* Y out of range */
      fabs(X)/d>(N+1)/2 ||                        /* X out of range */
      ((N%2)==0 && (xp<0.5-ap || xp>0.5+ap)) ||   /* Even N slit passage */
      ((N%2)==1 && (xp>ap && xp<1-ap)) ||         /* Odd  N slit passage */
      ((N%2)==1 && fabs(X)/d>N/2+0.5)             /* Fix at end for odd N */
    ) gptremoveparticle(par) ;

  return( 1 ) ;
}
