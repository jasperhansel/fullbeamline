/* sc3d.c - Space-Charge routine */

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

void spacecharge3D_init(gptinit *init)
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

/* Sim routine just forwards to sc3Dp2p */
static void sc_sim( double t, double *x, double *p, void *vinfo, int nOff, int nCPU )
{
  struct sc_info *info = (struct sc_info *)vinfo ;
  sc3Dp2p(t,x,p,info->zmin,info->zmax,nOff,nCPU) ;
}

/* Relativistic point-to-point interaction, also used by other elements as 'fallback' routine */
#define ffac (1/(4*gpt_pi*gpt_eps0))
void sc3Dp2p( double t, double *x, double *p, double zmin, double zmax, int nOff, int nCPU )
{
  struct par *pari, *parj ;

  double E[3] ;
  double r[3] ;

  for(int i=nOff ; i<numpar ; i+=nCPU)
    if( pars[i].alive && pars[i].Wr[2] > zmin && pars[i].Wr[2] < zmax )
      for(int j=0 ; j<numpar ; j++)
  {
    if( i==j || !pars[j].alive ) continue ;
    if( pars[j].Wr[2] < zmin || pars[j].Wr[2] > zmax ) continue ;

    pari = &pars[i] ;
    parj = &pars[j] ;

    double G = parj->G ;

    r[0] = pari->Wr[0] - parj->Wr[0] ;
    r[1] = pari->Wr[1] - parj->Wr[1] ;
    r[2] = pari->Wr[2] - parj->Wr[2] ;

    double cdot = INP( r, parj->GBr ) / (G+1.0) ;

    r[0] += cdot*parj->GBr[0] ;
    r[1] += cdot*parj->GBr[1] ;
    r[2] += cdot*parj->GBr[2] ;
	double r2 = VECSQR(r) + parj->r2 ;
	double r3 = r2 * sqrt(r2) ;

    double mfac = ffac*parj->q*parj->n / r3 ;
    E[0] = mfac*r[0] ;
    E[1] = mfac*r[1] ;
    E[2] = mfac*r[2] ;

    cdot = INP( E, parj->GBr ) / (G+1.0) ;

    pari->WE[0] += G*E[0] - cdot*parj->GBr[0] ;
    pari->WE[1] += G*E[1] - cdot*parj->GBr[1] ;
    pari->WE[2] += G*E[2] - cdot*parj->GBr[2] ;

    pari->WB[0] += (parj->GBr[1]*E[2] - parj->GBr[2]*E[1]) / gpt_c ;
    pari->WB[1] += (parj->GBr[2]*E[0] - parj->GBr[0]*E[2]) / gpt_c ;
    pari->WB[2] += (parj->GBr[0]*E[1] - parj->GBr[1]*E[0]) / gpt_c ;
  }
}
