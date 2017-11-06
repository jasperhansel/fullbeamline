/* coppscat.c - Copper scatter model */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "elem.h"

#define EMIN  50.0    /* in eV */
#define ETOP 400.0    /* in eV */
#define E3 800.0      /* in eV */
#define E4 3000.0     /* in eV */
#define fTOP 0.95
#define f3 0.45
#define f4 0.05
#define EREDUC 0.7      /* Energy reduction factor */
#define ESPREAD 0.5      /* Random energy spread */

#define EXTRAPOLATE 1e-6


struct copper_info
{
  double nmin, nmax ;

  struct scatter_info scatinfo ;
} ;

static void copper_scat(gptpar *par,double t,double dt, gpttrajectory *ptraj,void *info) ;
static double energydepfunc(double Ein) ;


void copperscatter_init(gptinit *init)
{
  struct copper_info *info ;
  char *name ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=3 )
    gpterror( "Syntax: %s(ECS,name,nmin,nmax)\n", gptgetname(init) ) ;

  info = (struct copper_info *)gptmalloc( sizeof(struct copper_info) ) ;

  name       = gptgetargstring(init,1) ;
  info->nmin = gptgetargdouble(init,2) ;
  info->nmax = gptgetargdouble(init,3) ;

  /* Install all functions */
  gptscatterinit(init,&info->scatinfo,name) ;
  gptinstallscatterfnc(name,copper_scat,info) ;
}


static void copper_scat(gptpar *par,double t,double dt, gpttrajectory *ptraj,void *vinfo)
{
  struct copper_info *info = (struct copper_info *)vinfo ;

  double forwardn, backwardn, Ein, energydep, angledep, Efac, ereduc ;
  double Gnew, lenGBnew, extratime ;
  double newdr[3], GBnew[3], vnew[3], rnew[3] ;
  int i ;

  /* Set output buffers and remove particles below absorption energy */
  gptscatterinitpar(&info->scatinfo) ;
  Ein = fabs(par->m*gpt_c*gpt_c/par->q) * (ptraj->Gint-1.0) ;
  if( Ein<=EMIN ) goto remove ;

  /* Calculate reflection chance */
  energydep = energydepfunc(Ein) ;
  angledep = 3-2*fabs(ptraj->inp) ;

  /* Scatter single particle with random generator */
  forwardn = backwardn = 0.0 ;
  if( rand()%3==0 )
    forwardn  = angledep*energydep*par->n ;
  else
    backwardn = angledep*energydep*par->n ;

  /* Scatter two particles */
  /* forwardn  = (1.0/3.0) * angledep*energydep*par->n ; */
  /* backwardn = (2.0/3.0) * angledep*energydep*par->n ; */

  /* Don't generate too small or too large particles */
  if( forwardn  < info->nmin ) forwardn = 0.0 ;
  if( backwardn < info->nmin ) backwardn = 0.0 ;
  if( forwardn  > info->nmax ) forwardn  = info->nmax ;
  if( backwardn > info->nmax ) backwardn = info->nmax ;
  if( forwardn==0.0 && backwardn==0.0 ) goto remove ;

  /* Calculate scattered particle parameters */
  Efac = par->n/(forwardn+backwardn) ;
  if( Efac > 1 ) Efac = 1.0 ;    

  ereduc = EREDUC+(ESPREAD*(double)rand()/RAND_MAX-ESPREAD/2) ;
  Gnew = 1.0+Efac*ereduc*(ptraj->Gint-1.0) ; /* Particles start with reduced energy */
  lenGBnew = sqrt(Gnew*Gnew-1.0) ;
/*  extratime = (1-trajectory->lambda+EXTRAPOLATE)*dt ; */
  extratime = EXTRAPOLATE*dt ;
  if( extratime<0 ) extratime=0 ;

  /* Forward scatter */
  if( forwardn!=0.0 )
  {
    for(i=0 ; i<3 ; i++) newdr[i] = ptraj->ndr[i]-2*ptraj->inp*ptraj->n[i] ;
    for(i=0 ; i<3 ; i++) GBnew[i] = newdr[i] * lenGBnew ;
    for(i=0 ; i<3 ; i++) vnew[i] = gpt_c*GBnew[i]/Gnew ;
    for(i=0 ; i<3 ; i++) rnew[i] = ptraj->P[i]+vnew[i]*extratime ;

    gptscatteraddparmqn(&info->scatinfo, gptgetparset("forwardscatter"), rnew, GBnew, par->m, par->q, forwardn ) ;
  }

  /* Backscatter */
  if( backwardn!=0.0 )
  {
    for(i=0 ; i<3 ; i++) newdr[i] = -ptraj->ndr[i] ;
    for(i=0 ; i<3 ; i++) GBnew[i] = newdr[i] * lenGBnew ;
    for(i=0 ; i<3 ; i++) vnew[i] = gpt_c*GBnew[i]/Gnew ;
    for(i=0 ; i<3 ; i++) rnew[i] = ptraj->P[i]+vnew[i]*extratime ;

    gptscatteraddparmqn(&info->scatinfo, gptgetparset("forwardscatter"), rnew, GBnew, par->m, par->q, backwardn ) ;
  }

  /* Remove original particle and exit */
remove:
  gptscatterremoveparticle(&info->scatinfo,ptraj,par) ;

  return ;
}


static double energydepfunc(double Ein)
{
  if( Ein<EMIN ) return(0.0 ) ;
  else if( Ein<ETOP ) return( fTOP/(ETOP-EMIN)*(Ein-EMIN) ) ;
  else if( Ein<E3 ) return( (f3-fTOP)/(E3-ETOP)*(Ein-ETOP) + fTOP ) ;
  else if( Ein<E4 ) return( (f4-f3)/(E4-E3)*(Ein-E3) + f3 ) ;
  return(f4) ;
}
