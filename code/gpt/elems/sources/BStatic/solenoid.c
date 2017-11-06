/* solenoid.c - Solenoid */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

struct solenoid_info
{
  double mu0I ;
  double R ;
} ;

static int solenoid_sim(gptpar *par,double t,struct solenoid_info *info) ;

static double rf(double x, double y, double z ) ;
static double rd(double x, double y, double z ) ;
static double vecxy( double x, double y ) ;


void solenoid_init(gptinit *init)
{
  struct solenoid_info *info ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=2 )
    gpterror( "Syntax: %s(ECS,R,I)\n", gptgetname(init) ) ;

  info = (struct solenoid_info *)gptmalloc( sizeof(struct solenoid_info) ) ;

  info->R = gptgetargdouble(init,1) ;
  info->mu0I = gpt_mu0*gptgetargdouble(init,2) ;

  gptaddEBelement( init, solenoid_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}


static int solenoid_sim(gptpar *par,double t,struct solenoid_info *info)
{
  double ra, tx, ty, tz ;
  double p, s2, rootp, sqrtrootp, rootm,  kp, rdv, rfv  ;
  double Bt ;

  ra=1/info->R ;
  tx=ra*X ;
  ty=ra*Y ;
  tz=ra*Z ;
  p=2*sqrt(tx*tx+ty*ty) ;
  s2=tx*tx+ty*ty+tz*tz ;
  rootp=1+p+s2 ;
  sqrtrootp=sqrt(rootp) ;
  rootm=1-p+s2 ;
  kp=rootm/rootp ;
  rdv=rd(0,kp,1) ;
  rfv=rf(0,kp,1) ;

  Bt=(tz*ra/(gpt_pi))*(rfv-2.0*rdv*(1+s2)/(3.0*rootp))/(rootm*sqrtrootp) ;
  BX=info->mu0I*vecxy(tx,ty)*Bt ;
  BY=info->mu0I*vecxy(ty,tx)*Bt ;
  BZ=info->mu0I*(ra/(2*gpt_pi))*((2-p)*rfv-2.0*rdv*p*(1-s2)/(3.0*rootp))/(rootm*sqrtrootp) ;

  return( 1 ) ;
}


static double vecxy( double x, double y )
{
  double tmp ;

  if( x==0 ) return( 0 ) ;

  if( fabs(x)>fabs(y) )
  {
    tmp=y/x ;
    return( (x>0 ? 1.0 : -1.0)/sqrt(1+tmp*tmp) ) ;
  } else
  {
    tmp=x/y ;
    return( (y>0 ? 1.0 : -1.0)*tmp/sqrt(1+tmp*tmp) ) ;
  }
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
