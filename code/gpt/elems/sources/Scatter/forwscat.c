/* forwscat.c - Forward scatter particles on boundary elements */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "elem.h"

#define EXTRAPOLATE 1e-6


struct forwscat_info
{
  double P ;
  double nmin ;

  struct scatter_info scatinfo ;
} ;

void forwscat_scat(gptpar *par,double t,double dt, gpttrajectory *trajectory,void *info) ;


void forwardscatter_init(gptinit *init)
{
  struct forwscat_info *info ;
  char *name ;
  int numarg ;

  gptbuildECS( init ) ;

  numarg = gptgetargnum(init) ;

  if( numarg!=2 && numarg!=3 )
    gpterror( "Syntax: %s(ECS,name,P,[nmin])\n", gptgetname(init) ) ;

  info = (struct forwscat_info *)gptmalloc( sizeof(struct forwscat_info) ) ;

  name       = gptgetargstring(init,1) ;
  info->P    = gptgetargdouble(init,2) ;

  if( numarg==3 )
    info->nmin = gptgetargdouble(init,3) ;
  else
    info->nmin = 0.0 ;

  if( info->P<0 || info->P>1 )
    gpterror( "Probability P must be between 0 and 1." ) ;
  if( info->nmin < 0 )
    gpterror( "Nmin must be nonnegative." ) ;

  /* Install all functions */
  gptscatterinit(init,&info->scatinfo,name) ;
  gptinstallscatterfnc(name,forwscat_scat,info) ;
}


void forwscat_scat(gptpar *par,double t,double dt, gpttrajectory *ptraj,void *vinfo)
{
  struct forwscat_info *info = (struct forwscat_info *)vinfo ;

  double Gnew, lenGBnew, extratime ;
  double newdr[3], GBnew[3], vnew[3], rnew[3] ;
  int i ;

  /* Increase output buffer, set Qout and Eout to 0 */
  gptscatterinitpar(&info->scatinfo) ;
  if( info->P==0 ) goto remove_and_exit ;

  /* Generate new particles */
  Gnew = ptraj->Gint ;
  lenGBnew = gptVECLEN(ptraj->GBint) ;
/*  extratime = (1-ptray->lambda+EXTRAPOLATE)*dt ; */
  extratime = EXTRAPOLATE*dt ;
  if( extratime<0 ) extratime=0 ;

  /* Forward scatter */
  for(i=0 ; i<3 ; i++) newdr[i] = ptraj->ndr[i]-2*ptraj->inp*ptraj->n[i] ;
  for(i=0 ; i<3 ; i++) GBnew[i] = newdr[i] * lenGBnew ;
  for(i=0 ; i<3 ; i++) vnew[i] = gpt_c*GBnew[i]/Gnew ;
  for(i=0 ; i<3 ; i++) rnew[i] = ptraj->P[i]+vnew[i]*extratime ;

  gptscatteraddparmqn(&info->scatinfo, gptgetparset("forwardscatter"), rnew, GBnew, par->m, par->q, par->n*info->P ) ;

  /* Remove original particle and exit */
remove_and_exit:
  gptscatterremoveparticle(&info->scatinfo,ptraj,par) ;

  return ;
}
