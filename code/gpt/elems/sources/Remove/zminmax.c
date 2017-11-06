/* zminmax.c: Remove all particles before zmin or after zmax */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include "elem.h"

/* Info structure containing all relevant parameters for this element */
struct zminmax_info
{
  double Zmin ;
  double Zmax ;
} ;

/* Forward declaration of the routine calculating the electromagnetic fields */
static int zminmax_sim(gptpar *par,double t,struct zminmax_info *info) ;

/* Initialization routine */
void zminmax_init(gptinit *init)
{
  struct zminmax_info *info ;

  /* Read Element Coordinate System (ECS) from parameter list */
  gptbuildECS( init ) ;

  /* Print usage line when the number of parameters is incorrect */
  if( gptgetargnum(init)!=2 )
    gpterror( "Syntax: %s(ECS,Zmin,Zmax)\n", gptgetname(init) ) ;

  /* Allocate memory for info structure */
  info = (struct zminmax_info *)gptmalloc( sizeof(struct zminmax_info) ) ;

  /* Read all parameters as doubles and store them in info structure */
  info->Zmin     = gptgetargdouble(init,1) ;
  info->Zmax     = gptgetargdouble(init,2) ;

  /* Register the routine calculating the electromagnetic fields to the GPT kernel */
  gptaddEBelement( init, zminmax_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}


/* The following routine calculates the electromagnetic fields */
static int zminmax_sim(gptpar *par,double t,struct zminmax_info *info)
{
  if( Z>=info->Zmin && Z<=info->Zmax ) return(1) ;

  gptremoveparticle(par) ;
  return(1) ;
}
