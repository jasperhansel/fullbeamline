/* map1D_TM.c - Reads a 1D table of axial Ez samples and a frequency and
 *              extrapolates these to a cylindrical symmetric field map
 *              for a TM mode cavity
 */

/* This file is part of GPT and belongs to:
 * S.B. van der Geer, M.J. de Loos and Pulsar Physics
 */

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include "elem.h"


struct map1d_info
{
  /* GDF information */
  struct gdfmem gm ;

  /* Coordinate and field information */
  double *z, *fz, *fz2 ;
  double omega, phi ;
  int points ;
} ;


static int map1d_sim(gptpar *par, double t, struct map1d_info *info) ;
static void map1d_exit(struct map1d_info *info) ;


void map1D_TM_init(gptinit *init)
{
  struct map1d_info *info ;
  struct gdfdata *ds ;

  /* Commandline parameters */
  char *zname ;
  char *fzname ;
  double ffac ;

  /* Grid characteristics */
  int i,points ;

  /* GDF information */
  struct gdff ingdff ;
  char *filename ;

  gptbuildECS( init ) ;

  /* Commandline parsing */
  if( gptgetargnum(init)!=6 )
    gpterror( "Syntax: %s(ECS,mapfile.gdf,z,Ez,fac,phi,omega)\n", gptgetname(init) ) ;

  info = (struct map1d_info *)gptmalloc( sizeof(struct map1d_info) ) ;

  /* Initialize GDF */
  filename = gptgetargstring(init,1) ;
  gdfsrinit( filename, &ingdff ) ;
  gdfrmem( &ingdff, &info->gm, GDFR_ABORTONERROR | GDFR_READARRAYS ) ;

  /* Retrieve z array */
  zname = gptgetargstring(init,2) ;
  ds=gptinputgetgroup(&info->gm,filename,zname) ;
  info->z=gptinputdoublearray(ds,zname,&points) ;
  info->points = points ;

  /* Retrieve field from datafile */
  fzname   = gptgetargstring(init,3) ;
  info->fz = gptinputdoublearraypoints(ds,fzname,points) ;

  /* Multiply field with extra factor */
  ffac = gptgetargdouble(init,4) ;
  for(i=0 ; i<points ; i++) info->fz[i] *= ffac ;

  /* Get frequency and phase */
  info->phi   = gptgetargdouble(init,5) ;
  info->omega = gptgetargdouble(init,6) ;

  /* Check ascending z-order */
  for( i=0 ; i<(points-1) ; i++ )
    if( info->z[i]>=info->z[i+1] )
      gpterror( "%s: z array must be in ascending order\n", filename ) ;

  /* Calculate spline information */
  info->fz2 = (double *)gptmalloc(points*sizeof(double)) ;
  if( spline( info->z, info->fz, info->fz2, points, NATSPLINE, NATSPLINE ) )
    gpterror( "%s: Error calculating spline\n", filename ) ;

  gptaddEBelement( init, map1d_sim, map1d_exit, GPTELEM_GLOBAL, info ) ;
}


static int map1d_sim(gptpar *par, double t, struct map1d_info *info)
{
  double fz, dfz ;
  double omega ;
  double phase, cosphase, sinphase ;
  double Ez, Eror, Bphior ;

  if( splintd( info->z, info->fz, info->fz2, info->points, Z, &fz, &dfz ) )
    return( 0 ) ;

  omega = info->omega ;
  phase = omega*t+info->phi ;
  cosphase = cos(phase) ;
  sinphase = sin(phase) ;

  Ez     =  cosphase*fz ;
  Eror   = -cosphase*dfz*0.5 ;
  Bphior = -sinphase*fz*omega/(2.0*gpt_c*gpt_c) ;

  EX =  X*Eror ;
  EY =  Y*Eror ;
  EZ =  Ez ;
  BX = -Y*Bphior ;
  BY =  X*Bphior ;

  return( 1 ) ;
}

static void map1d_exit(struct map1d_info *info)
{
  free(info->fz2) ;
  gdfmfreechilds( &info->gm.ds ) ;
  free(info) ;
}
