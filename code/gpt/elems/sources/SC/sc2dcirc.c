/* sc2dcirc.c - 2D Space-Charge routine using circle charges */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "elem.h"

struct sc_info
{
  double zmin ;
  double zmax ;
} ;

static void sc_sim( double t, double *x, double *p, void *info, int nOff, int nCPU ) ;
static void circlecharge( double *r, double R, double lambda, double *E ) ;
static double rf(double x, double y, double z ) ;
static double rd(double x, double y, double z ) ;

void spacecharge2Dcircle_init(gptinit *init)
{
  struct sc_info *info ;
  int numarg ;

  numarg = gptgetargnum(init) ;

  if( numarg!=0 && numarg!=2 )
    gpterror( "Syntax: %s([zmin,zmax])", gptgetname(init) ) ;

  info = (struct sc_info *)gptmalloc( sizeof(struct sc_info) ) ;

  if( numarg==2 )
  {
    info->zmin=gptgetargdouble(init,1) ;
    info->zmax=gptgetargdouble(init,2) ;
  } else
  {
    info->zmin=-DBL_MAX ;
    info->zmax= DBL_MAX ;
  }

  odemtaddfprfunction( ODEFNC_USR, sc_sim,info ) ;
}


static void sc_sim( double t, double *x, double *p, void *vinfo, int nOff, int nCPU )
{
  struct par *pari, *parj ;
  struct sc_info *info = (struct sc_info *)vinfo ;

  double E[3] ;
  double r[3] ;
  double Bcirc, Gcirc, Rcirc, lambda ;
  double zmin, zmax ;

  zmin    = info->zmin ;
  zmax    = info->zmax ;

  for(int i=nOff ; i<numpar ; i+=nCPU)
    if( pars[i].alive && pars[i].Wr[2] > zmin && pars[i].Wr[2] < zmax )
      for(int j=0 ; j<numpar ; j++)
  {
    if( i==j || !pars[j].alive ) continue ;
    if( pars[j].Wr[2] < zmin || pars[j].Wr[2] > zmax ) continue ;

    pari = &pars[i] ;
    parj = &pars[j] ;

    Bcirc = parj->GBr[2]/parj->G ;
    Gcirc = 1/sqrt(1-Bcirc*Bcirc) ;
    Rcirc = sqrt(parj->Wr[0]*parj->Wr[0]+parj->Wr[1]*parj->Wr[1]) ;

    r[0] = pari->Wr[0] ;
    r[1] = pari->Wr[1] ;
    r[2] = Gcirc*(pari->Wr[2] - parj->Wr[2]) ;

    lambda = parj->q*parj->n/(2*gpt_pi*Rcirc) ;
    circlecharge(r,Rcirc,lambda,E) ;

    pari->WE[0] += Gcirc*E[0] ;
    pari->WE[1] += Gcirc*E[1] ;
    pari->WE[2] += E[2] ;

    pari->WB[0] +=-Gcirc*Bcirc*E[1] / gpt_c ;
    pari->WB[1] += Gcirc*Bcirc*E[0] / gpt_c ;
    pari->WB[2] += 0.0 ;
  }
}

static void circlecharge( double *r, double R, double lambda, double *E )
{
  double rxy, Z2, rmR, d2, tmp2, tmp, ap, lamoeps0 ;
  double rfv, rdvo3 ;
  double Er ;

  rxy = sqrt(r[0]*r[0]+r[1]*r[1]) ;
  Z2 = r[2]*r[2] ;
  rmR = rxy-R ;
  d2 = rmR*rmR + Z2 ;
  if( d2==0 ) return ;
  tmp2 = d2+4*rxy*R ;
  tmp = sqrt(tmp2) ;
  ap = d2/tmp2 ;
  rfv=rf(0,ap,1) ;
  rdvo3=rd(0,ap,1)/3 ;
  lamoeps0=lambda/gpt_eps0 ;

  Er = lamoeps0*R/(gpt_pi*tmp)*(rmR*rfv/d2 + 2*R*(R*R-rxy*rxy+Z2)*rdvo3/(tmp2*d2)) ;
  gptr2carth(Er,r[0],r[1],&E[0],&E[1]) ;
  E[2] = lamoeps0*R*r[2]*(rfv-(1-ap)*rdvo3)/(gpt_pi*d2*tmp) ;
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
