/* iris.c - Kill particle on a foil with circular opening and generate new reflected particles */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "elem.h"

#define EXTRAPOLATE 1e-6

struct iris_info
{
  double Rin2, Rout2 ;
} ;

static int iris_bound(gptpar *par, double t, double dt, gpttrajectory *trajectory,struct iris_info *info ) ;

void scatteriris_init(gptinit *init)
{
  struct iris_info *info ;
  double R, tmp ;

  gptbuildECS( init ) ;

  if( gptgetargnum(init)!=2 )
    gpterror( "Syntax: %s(ECS,Rin,Rout)\n", gptgetname(init) ) ;

  info = (struct iris_info *)gptmalloc( sizeof(struct iris_info) ) ;

  R = gptgetargdouble(init,1) ; info->Rin2  = R*R ;
  R = gptgetargdouble(init,2) ; info->Rout2 = R*R ;

  /* Swap Rin and Rout when in reverse order */
  if( info->Rin2 > info->Rout2 )
  {
    tmp = info->Rout2 ;
    info->Rout2 = info->Rin2 ;
    info->Rin2 = tmp ;
  }

  gptaddboundaryelement( init, iris_bound, gptfree, 0, info ) ;
}


static int iris_bound(gptpar *par, double t, double dt, gpttrajectory *trajectory, struct iris_info *info )
{
  int i ;
  double *dr ;
  double lambda, PR2, Rstart2, Rend2 ;
  double P[3] ;

  /* See if the particle makes any chance of crossing the iris */
  if( trajectory->rstart[2] < 0.0 && trajectory->rend[2] < 0.0 ) return( 0 ) ;
  if( trajectory->rstart[2] > 0.0 && trajectory->rend[2] > 0.0 ) return( 0 ) ;

  Rstart2= trajectory->rstart[0]*trajectory->rstart[0]+trajectory->rstart[1]*trajectory->rstart[1] ;
  Rend2  = trajectory->rend[0]*trajectory->rend[0]+trajectory->rend[1]*trajectory->rend[1] ;
  if( Rstart2 < info->Rin2 && Rend2 < info->Rin2 ) return( 0 ) ;

  /* Calculate trajectory direction dr=rend-rstart */
  dr = trajectory->dr ;

  /* Calculate intersection point (P=rstart+lambda*dr, see Mathematica document */
  lambda = -trajectory->rstart[2]/dr[2] ;
  if( lambda > trajectory->lambda ) return( 0 ) ;
  for(i=0 ; i<3 ; i++) P[i] = trajectory->rstart[i] + lambda*dr[i] ;

  /* Test if P is within iris dimensions */
  PR2 = P[0]*P[0] + P[1]*P[1] ;
  if( PR2<info->Rin2 || PR2>info->Rout2 ) return( 0 ) ;

  /* Store result */
  for(i=0 ; i<3 ; i++) trajectory->P[i] = P[i] ;
  trajectory->lambda = lambda ;
  trajectory->n[0] =  0.0 ;
  trajectory->n[1] =  0.0 ;
  trajectory->n[2] = -1.0 ;  

  return( 1 ) ;
}
