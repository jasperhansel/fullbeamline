/* circlech.c - circlecharge */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

struct circlecharge_info
{
  double lamoeps0 ;
  double R ;
} ;

static int circlecharge_sim(gptpar *par,double t,struct circlecharge_info *info) ;
static double rf(double x, double y, double z ) ;
static double rd(double x, double y, double z ) ;


void circlecharge_init(gptinit *init)
{
  struct circlecharge_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=2 )
    gpterror( "Syntax: %s(ECS,R,lambda)\n", gptgetname(init) ) ;

  info = (struct circlecharge_info *)gptmalloc( sizeof(struct circlecharge_info) ) ;

  info->R = gptgetargdouble(init,1) ;
  info->lamoeps0 = gptgetargdouble(init,2)/gpt_eps0 ;

  if( info->R <= 0)
    gpterror( "%s: Radius must be positive\n", gptgetname(init) ) ;

  gptaddEBelement( init, circlecharge_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}


static int circlecharge_sim(gptpar *par,double t,struct circlecharge_info *info)
{
  double R, lamoeps0, r, Z2, rmR, d2, tmp2, tmp, ap ;
  double rfv, rdvo3 ;
  double Er ;

  R = info->R ;
  lamoeps0 = info->lamoeps0 ;
  r = sqrt(X*X+Y*Y) ;
  Z2 = Z*Z ;
  rmR = r-R ;
  d2 = rmR*rmR + Z2 ;
  if( d2==0 ) return( 0 ) ;
  tmp2 = d2+4*r*R ;
  tmp = sqrt(tmp2) ;
  ap = d2/tmp2 ;
  rfv=rf(0,ap,1) ;
  rdvo3=rd(0,ap,1)/3 ;

  Er = lamoeps0*R/(gpt_pi*tmp)*(rmR*rfv/d2 + 2*R*(R*R-r*r+Z2)*rdvo3/(tmp2*d2)) ;
  gptr2carth(Er,X,Y,&EX,&EY) ;
  EZ = lamoeps0*R*Z*(rfv-(1-ap)*rdvo3)/(gpt_pi*d2*tmp) ;

  return( 1 ) ;
}


/* From: Numerical recipes */
#define ERRRF 0.0020	/* MODIFIED for better accuracy */
static double rf(double x, double y, double z )
{
  double alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt ;

  xt=x ;
  yt=y ;
  zt=z ;
  do
  {
    sqrtx=sqrt(xt) ;
    sqrty=sqrt(yt) ;
    sqrtz=sqrt(zt) ;
    alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz ;
    xt=0.25*(xt+alamb) ;
    yt=0.25*(yt+alamb) ;
    zt=0.25*(zt+alamb) ;
    ave=(xt+yt+zt)/3.0 ;
    delx=(ave-xt)/ave ;
    dely=(ave-yt)/ave ;
    delz=(ave-zt)/ave ;
  } while( fabs(delx)>ERRRF || fabs(dely)>ERRRF || fabs(delz)>ERRRF ) ;
  e2=delx*dely-delz*delz ;
  e3=delx*dely*delz ;
  return( 1.0+(e2/24.0-0.1-3.0*e3/44.0)*e2+e3/14.0)/sqrt(ave) ;
}

/* From: Numerical recipes */
#define ERRRD 0.0015
static double rd(double x, double y, double z )
{
  double alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt ;

  xt=x ;
  yt=y ;
  zt=z ;
  sum=0.0 ;
  fac=1.0 ;
  do
  {
    sqrtx=sqrt(xt) ;
    sqrty=sqrt(yt) ;
    sqrtz=sqrt(zt) ;
    alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz ;
    sum+=fac/(sqrtz*(zt+alamb)) ;
    fac=0.25*fac ;
    xt=0.25*(xt+alamb) ;
    yt=0.25*(yt+alamb) ;
    zt=0.25*(zt+alamb) ;
    ave=0.2*(xt+yt+3.0*zt) ;
    delx=(ave-xt)/ave ;
    dely=(ave-yt)/ave ;
    delz=(ave-zt)/ave ;
  } while( fabs(delx)>ERRRD || fabs(dely)>ERRRD || fabs(delz)>ERRRD ) ;
  ea=delx*dely ;
  eb=delz*delz ;
  ec=ea-eb ;
  ed=ea-6.0*eb ;
  ee=ed+ec+ec ;
  return( 3.0*sum+fac*(1.0+ed*(-3.0/14.0+9.0*ed/88.0-4.5*delz*ee/26.0)
    +delz*(ee/6.0+delz*(-9.0*ec/22.0+3.0*delz*ea/26.0)))/(ave*sqrt(ave)) ) ;
}
