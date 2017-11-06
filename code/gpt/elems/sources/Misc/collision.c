/* collision.c - Detect particle-particle collisions and start a new "combined" particle */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

struct collision_info
{
  gptparset *ppp ;
} ;

static void collision_end(double tstart, double tend, double *dt, double *xstart, double *xend, void *vinfo, void *stepinfo) ;

void collision_init(gptinit *init)
{
  struct collision_info *info ;

  /* Print usage line when the number of parameters is incorrect */
  if( gptgetargnum(init)!=0 )
    gpterror( "Syntax: %s()\n", gptgetname(init) ) ;

  info = (struct collision_info *)gptmalloc( sizeof(struct collision_info) ) ;
  info->ppp = gptgetparset("wcs") ;

  odeaddendfunction( ODEFNC_USR, collision_end, info ) ;
}

static void collision_end(double tstart, double tend, double *dt, double *xstart, double *xend, void *vinfo, void *stepinfo)
{
  struct collision_info *info = (struct collision_info *)vinfo ;
  unsigned int k ;
  double l2 ;
  double parinm, parjnm, sumnm, newWr[3], newWGBr[3], newm, newq, newn ;

  /* Loop over all alive particles that will not be removed */
  for(int i=0 ; i<numpar ; i++)
    if( pars[i].alive && !pars[i].tokill )
      for(int j=0 ; j<numpar ; j++)
  {
    if( i==j || !pars[j].alive || pars[j].tokill ) continue ;

    /* Calculate distance^2 between the two particles, manual unroll */
    l2 = (pars[i].Wr[0]-pars[j].Wr[0])*(pars[i].Wr[0]-pars[j].Wr[0]) +
         (pars[i].Wr[1]-pars[j].Wr[1])*(pars[i].Wr[1]-pars[j].Wr[1]) +
         (pars[i].Wr[2]-pars[j].Wr[2])*(pars[i].Wr[2]-pars[j].Wr[2]) ;

    /* Test collision */
    if( l2 < pars[j].r2 )
    {
      /* Collision: initialization */
      parinm = pars[i].n*pars[i].m ;
      parjnm = pars[j].n*pars[j].m ;
      sumnm = parinm+parjnm ;
      if( sumnm==0.0 ) terminate( "Collision: Error: You must specify nmacro\n" ) ;

      /* Calaculate mass-weigthed new position and momentum */
      for(k=0;k<3;k++) newWr[k]  = (parinm*pars[i].Wr[k] +parjnm*pars[j].Wr[k] )/sumnm ;
      for(k=0;k<3;k++) newWGBr[k]= (parinm*pars[i].GBr[k]+parjnm*pars[j].GBr[k])/sumnm ;
      
      /* Calculate new mqn */
      newn = (pars[i].n+pars[j].n) ;
      newq = (pars[i].n*pars[i].q + pars[j].n*pars[j].q)/newn ;
      newm = (pars[i].n*pars[i].m + pars[j].n*pars[j].m)/newn ; /* Could use parinm etc. */

      /* Remove both particles */
      gptremoveparticle((&pars[i])) ; /* Double brackets due to error in macro def. */
      gptremoveparticle((&pars[j])) ;

      /* Start new particle */
      gptaddparmqn(info->ppp,newWr,newWGBr,newm,newq,newn) ;
    }
  }
}
